#pragma once
#include <vector>
#include <functional>

struct cubicSpline {
    double a, b, c, d;
};

std::vector<cubicSpline> findSplineCoefficients(const std::vector<double> &x,
                                                const std::function<double(double)> &f,
                                                double mu1,
                                                double mu2);
