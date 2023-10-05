#include "worker.h"
#include <QThread>
using namespace std;

Worker::Worker(QObject* parent): QObject(parent) {}

solution_t Worker::solve_sle_chebyshev_custom_region(
            double a, double b, double c, double d,
            size_t n, size_t m,
            const std::function<double(double, double)> &mu,
            const std::function<double(double, double)> &F,
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


    vector<double> tau(K);  // массив чебышёвского набора параметров
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
    if ((K & (K-1)) == 0) { // степень двойки
        std::vector<int> order1(K), order2(K);
        order1[0] = 0;
        for (int r = 1; r < K; r *= 2) {
            for (int i = 0; i < r; ++i) {
                order2[2*i] = order1[i];
                order2[2*i + 1] = 2*r - 1 - order1[i];
            }
            std::swap(order1, order2);
        }
        std::vector<double> tau_reorder(K);
        for (int i = 0; i < K; ++i)
            tau_reorder[i] = tau[order1[i]];
        swap(tau_reorder, tau);
    }

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
    double norm = 0.0;  // норма разности решений на соседних итерациях

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
    double r0_2 = 0.0;

    // предварительный цикл для расчёта ||r^(0)||
#pragma omp parallel for schedule(dynamic) reduction(+:r0_2)
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
            r0_2 += (*rc) * (*rc);
            ++vc;
            ++fc;
            ++rc;
        }
    }
    r0_2 = sqrt(r0_2);

    // основной цикл
    for (it = 0; it < iterations; ++it) {
        for (size_t s = 0; s < K; ++s) {
#pragma omp parallel
            {
                double ctau = tau[s];

#pragma omp for schedule(dynamic, 128)
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
            emit updateProgress("Выполнено циклов: " + QString::number(it) + ", точность: " + QString::number(norm));
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

    soln.n = n;
    soln.m = m;
    soln.k = K;
    soln.vVec = vvec;
    soln.eps = norm;
    soln.eps_max = epsilon;
    soln.r_max = r_max;
    soln.r_2 = sqrt(r_2);
    soln.iteration_num = it;
    soln.iteration_num_max = iterations;
    soln.row_bound = full_row_bound;
    return soln;
}

double u(double x, double y) {
    double a = sin(M_PI*x*y);
    return exp(a*a);
}

double f(double x, double y) {
    double arg = 2.0*M_PI*x*y;
    double a = sin(M_PI*x*y);
    return -M_PI*M_PI*(x*x + y*y)*exp(a*a)*(sin(arg)*sin(arg) + 2.0*cos(arg));
}

void Worker::launchMethod(size_t n, size_t m, size_t K, double eps, size_t nmax) {
    double h = (b-a)/n, k = (d-c)/m;
    double sinmaxh = sin(M_PI*(n-1)/(2*n));
    double sinmaxk = sin(M_PI*(m-1)/(2*m));
    double Mmax = 4.0 * ((sinmaxh*sinmaxh)/(h*h) + (sinmaxk*sinmaxk)/(k*k)),
           Mmin = 15.0;

    solution_t soln = solve_sle_chebyshev_custom_region(a, b, c, d, n, m, u, f, nmax, eps, K, Mmax, Mmin);
    emit sendResult(soln);
}

