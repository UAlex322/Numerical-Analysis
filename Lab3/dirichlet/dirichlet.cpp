#include "dirichlet.h"
#include <iostream>
#include <algorithm>
#include <cmath>

using std::max;
using std::abs;
using std::sqrt;

vector<double> solve_dirichlet_optimized(double a,
                                         double b,
                                         double c,
                                         double d,
                                         size_t n,
                                         size_t m,
                                         const function<double(double, double)> &m1,
                                         const function<double(double, double)> &m2,
                                         const function<double(double, double)> &m3,
                                         const function<double(double, double)> &m4,
                                         const function<double(double, double)> &f,
                                         size_t iterations,
                                         int64_t &step_num,
                                         double eps_max,
                                         double &eps) {
    vector<double> res((n+1)*(m+1));
    double *v = res.data();
    double h = (b - a)/n;
    double k = (d - c)/m;
    double x;
    double y;
    size_t N = n+1;

    x = a;
    for (size_t i = 0; i <= n; ++i, x += h) {
        v[i]       = m3(x, c);
        v[N*m + i] = m4(x, d);
    }
    y = c;
    for (size_t j = 0; j <= m; ++j, y += k) {
        v[N*j]     = m1(a, y);
        v[N*j + n] = m2(b, y);
    }

    double *g = new double[(n+1)*(m+1)];
    double *gcurr, *vcurr;
    double c1 =  0.5*k*k/(h*h + k*k);
    double c2 =  0.5*h*h/(h*h + k*k);
    double c3 =  0.5*h*h*k*k/(h*h + k*k);
    double new_value;  // новое значение в узле сетки
    double norm;       // норма разности решений на соседних итерациях

    y = c + k;
    gcurr = g + N + 1;  // g[1][1]
    for (size_t j = 1; j < m; ++j, y += k) {
        x = a + h;
        for (size_t i = 1; i < n; ++i, x += h) {
            *gcurr = c3 * f(x,y);
            ++gcurr;  // g[j][i] -> g[j][i+1]
        }
        gcurr += 2;   // g[j][n] -> g[j+1][1]
    }

    step_num = 0;
    for (size_t t = 0; t < iterations; ++t) {
        norm = 0.0;
        gcurr = g + N + 1;  // g[1][1]
        vcurr = v + N + 1;  // v[1][1]
        for (size_t j = 1; j < m; ++j) {
            for (size_t i = 1; i < n; ++i) {
                new_value = c1 * (*(vcurr - 1) + *(vcurr + 1)) + c2 * (*(vcurr - N) + *(vcurr + N)) + *gcurr;
                norm = max(norm, abs(new_value - *vcurr));
                *vcurr = new_value;
                ++vcurr;  // v[j][i] --> v[j][i+1]
                ++gcurr;  // g[j][i] --> g[j][i+1]
            }
            vcurr += 2;  // v[j][n] --> v[j+1][1]
            gcurr += 2;  // g[j][n] --> g[j+1][1] 
        }
        ++step_num;
        eps = norm;
        if (norm < eps_max) break;
    }

    delete[] g;
    return res;
}

double solution_error(const function<double(double, double)> &u, const vector<double> &v_vector, size_t n, size_t m) {
    const double *v = v_vector.data();
    double error = 0.0;
    double h = 1.0/n;
    double k = 1.0/m;
    for (size_t j = 1; j < m; ++j)
        for (size_t i = 1; i < n; ++i) {
            error = max(error, abs(u(h*i, k*j) - v[j*(n+1) + i]));
        }

    return error;
}

double discrepancy(const function<double(double, double)> &f,
                   const vector<double> &v_vector,
                   size_t n,
                   size_t m) {
    size_t N = n+1;
    const double *v = v_vector.data();
    double x;
    double y;
    double h = 1.0/n;
    double k = 1.0/m;
    double h2 = (double)(n*n);
    double k2 = (double)(m*m);
    double a2 = -2.0*(h2 + k2);
    double R_norm = 0.0;
    double diff;

    y = k;
    for (size_t j = 1; j < m; ++j, y += k) {
        x = h;
        for (size_t i = 1; i < n; ++i, x += h) {
            diff = a2 * v[N*j + i] + h2 * (v[N*j + i-1] + v[N*j + i+1]) + k2 * (v[N*(j-1) + i] + v[N*(j+1) + i]) + f(x,y);
            R_norm += diff * diff;
        }
    }
    
    return sqrt(R_norm);
}

vector<double> solve_dirichlet(double a,
                               double b,
                               double c,
                               double d,
                               size_t n,
                               size_t m,
                               const function<double(double, double)> &m1,
                               const function<double(double, double)> &m2,
                               const function<double(double, double)> &m3,
                               const function<double(double, double)> &m4,
                               const function<double(double, double)> &f,
                               size_t iterations,
                               double epsilon) {
    vector<double> res((n+1)*(m+1));
    double *v = res.data();
    double h = (b - a)/n;
    double k = (d - c)/m;
    double x;
    double y;
    size_t N = n+1;

    x = a;
    for (size_t i = 0; i <= n; ++i, x += h) {
        v[i]           = m3(x, c);
        v[(n+1)*m + i] = m4(x, d);
    }
    y = c;
    for (size_t j = 0; j <= m; ++j, y += k) {
        v[(n+1)*j]     = m1(a, y);
        v[(n+1)*j + n] = m2(b, y);
    }

    double *g = new double[(n+1)*(m+1)];
    double c1 = 0.5*k*k/(h*h + k*k);
    double c2 = 0.5*h*h/(h*h + k*k);
    double c3 = 0.5*h*h*k*k/(h*h + k*k);
    double new_value;
    double norm;

    y = c + k;
    for (size_t j = 1; j < m; ++j, y += k) {
        x = a + h;
        for (size_t i = 1; i < n; ++i, x += h)
            g[N*j + i] = c3 * f(x,y);
    }

    for (size_t t = 0; t < iterations; ++t) {
        norm = 0.0;
        y = c + k;

        for (size_t j = 1; j < m; ++j, y += k) {
            x = a + h;
            for (size_t i = 1; i < n; ++i, x += h) {
                new_value = c1 * (v[N*j + i-1] + v[N*j + i+1]) + c2 * (v[N*(j-1) + i] + v[N*(j+1) + i]) + g[N*j + i];
                norm = max(norm, abs(new_value - v[N*j + i]));
                v[N*j + i] = new_value;
            }
        }
        
        if (norm < epsilon) break;
    }

    delete[] g;
    return res;
}