#include "spline.h"
using std::vector;
using std::function;

std::vector<double> tridiagonalSolver(size_t n,
                                      const std::vector<double> &A,
                                      const std::vector<double> &B,
                                      const std::vector<double> &C,
                                      const std::vector<double> &pi,
                                      double kappa1,
                                      double kappa2,
                                      double m1,
                                      double m2) {
    // коэффициенты метода прогонки
    std::vector<double> alpha(n+1), beta(n+1);
    // вектор решения
    std::vector<double> y(n+1);
    // начальные значения коэффициентов метода прогонки
    alpha[1] = kappa1, beta[1] = m1;

    for (size_t i = 1; i < n; ++i) {
        // знаменатель в коэффициентах метода прогонки
        double divide = 1/(C[i] - alpha[i]*A[i]);
        alpha[i+1] = B[i] * divide;
        beta[i+1]  = (pi[i] + beta[i]*A[i]) * divide;
    }

    y[n] = (m2 + beta[n]*kappa2)/(1 - alpha[n]*kappa2);
    for (size_t i = n; i >= 1; --i)
        y[i-1] = alpha[i]*y[i] + beta[i];

    return y;
}


vector<cubicSpline> findSplineCoefficients(const vector<double> &x,
                                           const function<double(double)> &f,
                                           double mu1,
                                           double mu2) {
    size_t n = x.size() - 1;
    vector<double> diag1(n),
                   diag2(n),
                   diag3(n),
                   pi(n);
    vector<cubicSpline> spline(n); // вектор коэффициентов сплайнов
    vector<double> df(n+1), // df[i] = f(x[i+1]) - f(x[i])
                    h(n),  //  h[i] = x[i+1] - x[i]
                    c;     //  c[i] = P_i''(x[i])

    // заполняем h[] и df[]
    for (size_t i = 0; i < n; ++i)
        h[i] = x[i+1] - x[i];
    df[0] = f(x[0]);
    for (size_t i = 0; i < n; ++i) {
        spline[i].a = df[i+1] = f(x[i+1]); // 'a' вычисляется здесь
        df[i] = df[i+1] - df[i];
    }
    // заполняем диагонали матрицы и правую часть СЛАУ, решение которой - c[]
    for (size_t i = 1; i < n; ++i) {
        diag1[i] = h[i-1];
        diag2[i] = -2.0*(h[i-1] + h[i]);
        diag3[i] = h[i];
        pi[i] = -6.0*(df[i]/h[i] - df[i-1]/h[i-1]);
    }

    // находим c[], по его значениям - b[] и d[]
    c = tridiagonalSolver(n, diag1, diag3, diag2, pi, 0.0, 0.0, mu1, mu2);
    for (size_t i = 0; i < n; ++i) {
        spline[i].b = df[i]/h[i] + (c[i]/6.0 + c[i+1]/3.0)*h[i];
        spline[i].c = 0.5 * c[i+1];
        spline[i].d = (c[i+1] - c[i])/(6.0 * h[i]);
    }

    return spline;
}
