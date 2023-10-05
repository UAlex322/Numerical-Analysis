#ifndef DIRICHLET_H
#define DIRICHLET_H

#include <vector>
#include <functional>
#include <QLabel>
#include "gui_updater.h"

struct solution_t {
    std::vector<size_t> row_bound;  // границы строк, входящих в сетку
    std::vector<double> vVec;       // вектор численного решения
    double r_max = 0.0;             // максимум-норма невязки на выходе
    double r0_max = 0.0;            // максимум-норма невязки в начале
    double r_2 = 0.0;               // евклидова норма невязки на выходе
    double eps;                     // точность метода на выходе
    size_t iteration_num;           // число итераций на выходе
};

inline void discrepancy(const double *x,
                        const double *b,
                        double *r,
                        size_t m,
                        size_t n,
                        double h2,
                        double k2,
                        double a2);

solution_t solve_sle_chebyshev(double a,
                               double b,
                               double c,
                               double d,
                               size_t n,
                               size_t m,
                               const std::function<double(double, double)> &mu,
                               const std::function<double(double, double)> &F,
                               size_t iterations,
                               double epsilon,
                               size_t k,
                               double Mmax,
                               double Mmin);

solution_t solve_sle_chebyshev_custom_region(
                    double a, double b, double c, double d,
                    size_t n, size_t m,
                    const std::function<double(double, double)> &mu,
                    const std::function<double(double, double)> &F,
                    size_t iterations, double epsilon,
                    size_t K, double Mmax, double Mmin);

#endif // DIRICHLET_H

