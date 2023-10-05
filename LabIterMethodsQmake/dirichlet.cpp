#define _USE_MATH_DEFINES
#include "dirichlet.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <cmath>
#include <QWidget>
#include <QLabel>
#include <omp.h>

using std::vector;
using std::function;
using std::max;
using std::abs;

class Worker: public QObject {
    Q_OBJECT

    QWidget *progressWidget = new QWidget;
    QLabel *progressLabel = new QLabel(progressWidget);
    progressWidget->setGeometry(300, 300, 300, 100);
    progressLabel->setGeometry(20, 20, 260, 40);
    progressLabel->setText("Here we go!");
};

inline void discrepancy(const double *x,
                        const double *f,
                        double *r,
                        size_t m,
                        size_t n,
                        double h2,
                        double k2,
                        double a2) {
    size_t N = n+1;
    const double *xc;
    const double *fc;
    double *rc;

#pragma omp parallel for private(xc, fc, rc) schedule(dynamic)
    for (size_t j = 1; j < m; ++j) {
        xc = x + j*N + 1;
        fc = f + j*N + 1;
        rc = r + j*N + 1;
        for (size_t i = 1; i < n; ++i) {
             *rc = h2 * (*(xc - 1) + *(xc + 1)) +
                      k2 * (*(xc - N) + *(xc + N)) +
                     (*fc) - a2 * (*xc);
            ++xc;
            ++fc;
            ++rc;
        }
    }
}

solution_t solve_sle_chebyshev(double a,
                               double b,
                               double c,
                               double d,
                               size_t n,
                               size_t m,
                               const function<double(double, double)> &mu,
                               const function<double(double, double)> &F,
                               size_t iterations,
                               double epsilon,
                               size_t K,
                               double Mmax,
                               double Mmin) {
    vector<double> vvec((n+1)*(m+1));
    vector<double>  aux((n+1)*(m+1));
    vector<double> fvec((n+1)*(m+1));
    vector<double> rvec((n+1)*(m+1));
    vector<double> vsvec((n+1)*(m+1));

    double *v  = vvec.data();
    double *v1 =  aux.data();
    double *f  = fvec.data();
    double *r  = rvec.data();
    double *vs = vsvec.data();
    double *tau = new double[K];

    double h = (b - a)/n;
    double k = (d - c)/m;
    double x;
    double y;
    size_t N = n+1;

    // fill data

    for (int s = 0; s < K; ++s)
        tau[s] = 2.0 / ((Mmax + Mmin) + (Mmax - Mmin)*cos(M_PI*(1 + 2*s)/(2*K)));

    x = a;
    for (size_t i = 0; i <= n; ++i, x += h) {
        v1[i] = vs[i] = v[i] = mu(x, c);
        v1[N*m + i] = vs[N*m + i] = v[N*m + i] = mu(x, d);
    }
    y = c;
    for (size_t j = 0; j <= m; ++j, y += k) {
        v1[N*j] = vs[N*j] = v[N*j] = mu(a, y);
        v1[N*j + n] = vs[N*j + n] = v[N*j + n] = mu(b, y);
    }
    
    double *fc,
           *vc,
           *vnc,
           *vsc,
           *rc;
    double h2 = 1.0/(h*h);
    double k2 = 1.0/(k*k);
    double a2 = 2.0*(h2 + k2);
    double norm;    // норма разности решений на соседних итерациях

#pragma omp parallel for private(fc, x, y)
    for (size_t j = 1; j < m; ++j) {
        fc = f + j*N + 1;
        y = c + k*j;
        x = a + h;
        for (size_t i = 1; i < n; ++i, x += h) {
            *fc = F(x,y);
            ++fc;
        }
        y += k;
    }

    solution_t soln;
    size_t it;
    for (it = 0; it < iterations; ++it) {
        for (size_t s = 0; s < K; ++s) {
            discrepancy(v, f, r, m, n, h2, k2, a2);

#pragma omp parallel 
            {
                double ctau = tau[s];
#pragma omp for private(rc, vc, vnc) schedule(dynamic)
                for (size_t j = 1; j < m; ++j) {
                    rc = r + j*N + 1;
                    vc = v + j*N + 1;
                    vnc = v1 + j*N + 1;
                    for (size_t i = 1; i < n; ++i) {
                        *vnc = *vc + ctau * (*rc);
                        ++vc;
                        ++vnc;
                        ++rc;
                    }
                }
            }
            std::swap(v, v1); // v^(s) <- v^(s+1)
        }
        
        norm = 0.0;
#pragma omp parallel for private(vc, vsc) schedule(dynamic) reduction(max:norm)
        for (size_t j = 1; j < m; ++j) {
            vc  =  v + j*N + 1;
            vsc = vs + j*N + 1;
            for (size_t i = 1; i < n; ++i) {
                norm = std::max(norm, std::abs(*vc - *vsc));
                *vsc = *vc;
                ++vc;
                ++vsc;
            }
        }
        if (norm < epsilon)
            break;
    }
    delete[] tau;
    soln.iteration_num = std::min(iterations, it+1);
    soln.r_max = norm;
    soln.vVec = std::move(vsvec);
    return soln;
}

inline void discrepancy_custom_region(const double *x,
                                      const double *f,
                                      double *r,
                                      const vector<size_t> &row_bound,
                                      size_t m,
                                      size_t n,
                                      double h2,
                                      double k2,
                                      double a2) {
    size_t N = n+1;
    size_t shift;
    const double *xc;
    const double *fc;
    double *rc;

    for (size_t j = 1; j < m; ++j) {
        shift = N*j + row_bound[2*j];
        xc = x + shift;
        fc = f + shift;
        rc = r + shift;
        for (size_t i = row_bound[2*j]; i < row_bound[2*j + 1]; ++i) {
            *rc = h2 * (*(xc - 1) + *(xc + 1)) +
                k2 * (*(xc - N) + *(xc + N)) +
                (*fc) - a2 * (*xc);
            ++xc;
            ++fc;
            ++rc;

        }
    }
}

solution_t solve_sle_chebyshev_custom_region(
        double a, double b, double c, double d,
        size_t n, size_t m,
        const function<double(double, double)> &mu,
        const function<double(double, double)> &F,
        size_t iterations, double epsilon,
        size_t K, double Mmax, double Mmin) {
    vector<double>   vvec((n+1)*(m+1));
    vector<double> auxvec((n+1)*(m+1));
    vector<double>   fvec((n+1)*(m+1));
    vector<double>   rvec((n+1)*(m+1));
    vector<double>  vsvec((n+1)*(m+1));
    double *v  =     vvec.data();  // текущее приближение
    double *v1 =   auxvec.data();  // следующее приближение
    double *f  =     fvec.data();  // значениямя f(x_i, y_j)
    double *r  =     rvec.data();  // невязка
    double *vs =    vsvec.data();  // предыдущее приближение для
                                   // проверки критерия остановки по точности

    QWidget *progressWidget = new QWidget;
    QLabel *progressLabel = new QLabel(progressWidget);
    progressWidget->setGeometry(300, 300, 300, 100);
    progressLabel->setGeometry(20, 20, 260, 40);
    progressLabel->setText("Here we go!");
    progressWidget->show();

    double *tau = new double[K];  // массив чебышёвского набора параметров
    vector<size_t> row_bound(2*(m+1));  // границы строк внутренних узлов сетки
    vector<size_t> full_row_bound(2*(m+1));  // границы строк узлов всей сетки

    // заполнение границ строк внутри сетки
    for (size_t i = 0; i <= 2*(m/4); i += 2) {
        row_bound[i] = n/4 + 1;
        row_bound[i+1] = 3*n/4;
    }
    for (size_t i = 2*(m/4) + 2; i < 6*(m/4); i += 2) {
        row_bound[i] = 1;
        row_bound[i+1] = n;
    }
    for (size_t i = 6*(m/4); i <= 2*m; i += 2) {
        row_bound[i] = n/4 + 1;
        row_bound[i+1] = n;
    }

    for (size_t i = 0; i < 2*(m/4); i += 2) {
        full_row_bound[i] = n/4;
        full_row_bound[i+1] = 3*n/4 + 1;
    }
    for (size_t i = 2*(m/4); i <= 6*(m/4); i += 2) {
        full_row_bound[i] = 0;
        full_row_bound[i+1] = n+1;
    }
    for (size_t i = 6*(m/4) + 2; i <= 2*m; i += 2) {
        full_row_bound[i] = n/4;
        full_row_bound[i+1] = n+1;
    }

    // вычисление K параметров итерационного метода
    for (size_t s = 0; s < K; ++s)
        tau[s] = 2.0 / (Mmax * (1.0 + cos((M_PI*(1 + 2*s))/(2*K))) +
                        Mmin * (1.0 - cos((M_PI*(1 + 2*s))/(2*K))));

    int m14 = m/4;
    int m34 = 3*m/4;
    int n14 = n/4;
    int n34 = 3*n/4;

    // вид области:
    //   ______            
    //   |    |            
    // --     |            
    // |      |            
    // --    --            
    //   |__|
    
    double h = (b - a)/n;
    double k = (d - c)/m;
    double x;
    double y;
    size_t N = n+1;

    // заполнение горизонтальных границ
    for (size_t i = n/4; i <= 3*n/4; ++i)
        v1[i] = vs[i] = v[i] = mu(a + i*h, c);
    for (size_t i = 0; i <= n/4; ++i)
        v1[N*m14 + i] = vs[N*m14 + i] = v[N*m14 + i] = mu(a + i*h, c + m14*k);
    for (size_t i = 3*n/4; i <= n; ++i)
        v1[N*m14 + i] = vs[N*m14 + i] = v[N*m14 + i] = mu(a + i*h, c + m14*k);
    for (size_t i = 0; i <= n/4; ++i)
        v1[N*m34 + i] = vs[N*m34 + i] = v[N*m34 + i] = mu(a + i*h, c + m34*k);
    for (size_t i = n/4; i <= n; ++i)
        v1[N*m + i] = vs[N*m + i] = v[N*m + i] = mu(a + i*h, d);

    // заполнение вертикальных границ
    for (size_t j = m/4; j <= 3*m/4; ++j)
        v1[N*j] = vs[N*j] = v[N*j] = mu(a, c + j*k);
    for (size_t j = 0; j <= m/4; ++j)
        v1[N*j + n14] = vs[N*j + n14] = v[N*j + n14] = mu(a + n14*h, c + j*k);
    for (size_t j = 3*m/4; j <= m; ++j)
        v1[N*j + n14] = vs[N*j + n14] = v[N*j + n14] = mu(a + n14*h, c + j*k);
    for (size_t j = 0; j <= m/4; ++j)
        v1[N*j + n34] = vs[N*j + n34] = v[N*j + n34] = mu(a + n34*h, c + j*k);
    for (size_t j = m/4; j <= m; ++j)
        v1[N*j + n] = vs[N*j + n] = v[N*j + n] = mu(b, c + j*k);

    double h2 = 1.0/(h*h);
    double k2 = 1.0/(k*k);
    double a2 = 2.0*(h2 + k2);
    double norm;  // норма разности решений на соседних итерациях

    // заполнение 'f' значениями f(x_i, y_j)
    y = c + k;
    for (size_t j = 1; j < m; ++j) {
        double *fc = f + j*N + row_bound[2*j];
        x = a + row_bound[2*j] * h;
        for (size_t i = row_bound[2*j]; i < row_bound[2*j + 1]; ++i) {
            *fc = F(x,y);
            ++fc;
            x += h;
        }
        y += k;
    }


    solution_t soln; // структура с данными о решении
    size_t it;       // номер итерации метода


    // предварительный цикл для расчёта ||r^(0)||
#pragma omp parallel for schedule(dynamic)
    for (size_t j = 1; j < m; ++j) {
        size_t shift = j*N + row_bound[2*j];
        double *vc = v + shift; // строка начального приближения
        double *fc = f + shift; // строка 'f'
        double *rc = r + shift; // строка невязки

        for (size_t i = row_bound[2*j]; i < row_bound[2*j + 1]; ++i) {
            // r[i][j] = 1/h^2 * (v[i-1][j] + v[i+1][j])
            //         + 1/k^2 * (v[i][j-1] + v[i][j+1])
            //         + f[i][j] - 2(1/h^2 + 1/k^2) * v[i][j]
            *rc = h2 * (*(vc - 1) + *(vc + 1)) +
                  k2 * (*(vc - N) + *(vc + N)) +
                  (*fc) - a2 * (*vc);
            if (abs(*rc) > soln.r0_max)
                soln.r0_max = *rc;
            ++vc;
            ++fc;
            ++rc;
        }
    }

    // основной цикл
    for (it = 0; it < iterations; ++it) {
        for (size_t s = 0; s < K; ++s) {
#pragma omp parallel
            {
                double ctau = tau[s];

#pragma omp for schedule(dynamic, 64)
                for (size_t j = 1; j < m; ++j) {
                    size_t shift = j*N + row_bound[2*j];
                    size_t row_begin = row_bound[2*j];
                    size_t row_end = row_bound[2*j + 1];
                    double *vc = v + shift;    // строка текущего приближения
                    double *vnc = v1 + shift;  // строка следующего приближения
                    double *fc = f + shift;    // строка 'f'
                    double *rc = r + shift;    // строка невязки

                    for (size_t i = row_begin; i < row_end; ++i) {
                        *rc = h2 * (*(vc - 1) + *(vc + 1)) +
                              k2 * (*(vc - N) + *(vc + N)) +
                              (*fc) - a2 * (*vc);
                        *vnc = *vc + ctau * (*rc);

                        ++vc;
                        ++vnc;
                        ++fc;
                        ++rc;
                    }
                }
            }
            std::swap(v, v1); // v^(s) <- v^(s+1), экономия времени и ресурсов
        }

        // проверка критерия остановки по точности
        norm = 0.0;
#pragma omp parallel for reduction(max:norm) schedule(dynamic)
        for (size_t j = 1; j < m; ++j) {
            size_t shift = j*N + row_bound[2*j];
            double  *vc =  v + shift;  // строка текущего приближения
            double *vsc = vs + shift;  // строка предыдущего приближения

            for (size_t i = row_bound[2*j]; i < row_bound[2*j + 1]; ++i) {
                norm = std::max(norm, std::abs(*vc - *vsc));
                *vsc = *vc;
                ++vc;
                ++vsc;
            }
        }
        if (norm < epsilon || it % 10 == 0) {
            progressLabel->setText(QString::number(it));
        }
        if (norm < epsilon)
            break;
    }

    // расчёт норм невязки по окончании работы метода
    double r_max = 0.0;
    double r_2 = 0.0;
#pragma omp parallel for reduction(max:r_max) \
        reduction(+:r_2) schedule(dynamic)
    for (size_t j = 1; j < m; ++j) {
        size_t shift = j*N + row_bound[2*j];
        double *rc =  r + shift;
        for (size_t i = row_bound[2*j]; i < row_bound[2*j + 1]; ++i) {
            r_max = std::max(r_max, std::abs(*rc));
            r_2 += (*rc) * (*rc);
            ++rc;
        }
    }

    delete[] tau;
    delete progressWidget;
    soln.vVec = vvec;
    soln.eps = norm;
    soln.r_max = r_max;
    soln.r_2 = sqrt(r_2);
    soln.iteration_num = it;
    soln.row_bound = full_row_bound;
    return soln;
}
