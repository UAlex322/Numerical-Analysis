#ifndef DIRICHLET_H
#define DIRICHLET_H

#include <vector>
#include <functional>
using std::vector;
using std::function;


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
                                         double &eps);

double solution_error(const function<double(double, double)> &u, 
                      const vector<double> &v_vector, 
                      size_t n, 
                      size_t m);

double discrepancy(const function<double(double, double)> &f, 
                   const vector<double> &v_vector, 
                   size_t n, 
                   size_t m);

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
                               double epsilon);

#endif // DIRICHLET_H

