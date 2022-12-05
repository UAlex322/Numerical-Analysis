#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include <array>
#include <limits>
#include <utility>

////////////////////////////////////////////////////////////////////////////////////


// Точка dim-мерного пространства
template <size_t dim>
using Point = std::array<double, dim>;

template <size_t dim>
Point<dim> operator+(const Point<dim> &a, const Point<dim> &b) {
    Point<dim> result(a);

    for (size_t i = 0; i < dim; ++i)
        result[i] += b[i];

    return result;
}

template <size_t dim>
Point<dim> operator*(double x, const Point<dim> &a) {
    Point<dim> result(a);

    for (size_t i = 0; i < dim; ++i)
        result[i] *= x;

    return result;
}

template <size_t dim>
Point<dim> operator*(const Point<dim> &a, double x) {
    Point<dim> result(a);

    for (int i = 0; i < dim; ++i)
        result[i] *= x;

    return result;
}

// Расстояние между точками (нужно для оценки локальной погрешности)
template <size_t dim>
double dist(const Point<dim> &a, const Point<dim> &b) {
    double result = 0.0;
    for (size_t i = 0; i < dim; ++i)
        result = std::max(result, std::abs(b[i] - a[i]));
    return result;
}

template <>
double dist(const Point<1> &a, const Point<1> &b) {
    return std::abs(b[0] - a[0]);
}

double dist(double a, double b) {
    return std::abs(b - a);
}


////////////////////////////////////////////////////////////////////////////////////


// Тип функции, обозначающей правую часть ДУ
template <typename value_t>
using DE_function = std::function<double(double, value_t)>;

// Обозначение набора функций, задающих правую часть ДУ
template <typename value_t, size_t dim>
using Diff_equation = std::array<DE_function<value_t>, dim>;

// Тип функции, реализующей шаг метода
template <typename value_t, size_t dim>
using Method_function = std::function<value_t(double, value_t, Diff_equation<value_t, dim>, double)>;

// Данные о полученном решении ДУ
template <typename value_t, size_t dim>
struct Entry {
    double x;
    value_t v;
    value_t v2;
    double s;
    double h;
    int c_plus;
    int c_minus;
};


// Численный метод
template <typename value_t, size_t dim>
struct Method {
    const Method_function<value_t, dim> next_point; // Шаг метода для получения следующей точки
    const size_t p;			          // Порядок метода
    const double pow2;                // 2^p
    const double s_mult;              // 1/(2^p - 1)
    const double s_star_mult;         // 2^p/(2^p - 1)

    Method(const Method_function<value_t, dim> &next_point, size_t p):
        next_point(next_point),
        p(p),
        pow2(pow(2.0, p)),
        s_mult(1/(pow2 - 1)),
        s_star_mult(1.0 + s_mult)
    {}
};


////////////////////////////////////////////////////////////////////////////////////


// Общий вид метода без регулировки шага
template <typename value_t, size_t dim>
std::vector<Entry<value_t, dim>> ivp_no_step_adjust (const Method<value_t, dim> &method,
                                                     const Diff_equation<value_t, dim> &F,
                                                     const double x0,
                                                     const value_t &u0,
                                                     const double h,
                                                     const double x_max,
                                                     const double x_eps)
{
    double x = x0;
    value_t v = u0;
    value_t v_half;
    double h_half = 0.5 * h;
    std::vector<Entry<value_t,dim>> soln;

    soln.push_back({x,v,v,0.0,h,0,0});
    while (x + h < x_max - x_eps) {
        v_half = method.next_point(x + h_half, method.next_point(x, v, F, h_half), F, h_half);
        v = method.next_point(x, v, F, h);
        x += h;

        soln.push_back({x, v, v_half, method.s_star_mult * dist(v_half,v), h, 0, 0});
    }

    return soln;
}


// Общий вид метода с регулировкой шага
template <typename value_t, size_t dim>
std::pair<std::vector<Entry<value_t, dim>>, int> ivp_step_adjust (const Method<value_t, dim> &method,
                                                  const Diff_equation<value_t, dim> &F,
                                                  const double x0,
                                                  const value_t &u0,
                                                  double h,
                                                  const double x_max,
                                                  const double x_eps,
                                                  const double v_eps,
                                                  int nmax = 10000)
{
    int iteration = 0;
    int plus_count = 0;
    int minus_count = 0;

    double x = x0;
    value_t v = u0, v_next, v_half;
    double h_half = 0.5*h;
    double s;       // контрольная величина;
    double eps_lower_bound = v_eps/pow(2.0, method.p + 1); // при |S| < eps/2^(p+1) увеличивается шаг

    std::vector<Entry<value_t, dim>> soln;
    soln.push_back({x0,u0,u0,0.0,h,0,0});

    while (iteration < nmax && x < x_max - x_eps) {
        v_half = method.next_point(x + h_half, method.next_point(x, v, F, h_half), F, h_half);
        v_next = method.next_point(x, v, F, h);
        s      = method.s_mult * dist(v_half, v_next);

        if (std::abs(s) > v_eps || x + h > x_max - x_eps) {
            h = h_half;
            h_half *= 0.5;
            ++minus_count;
        }
        else {
            x += h;
            if constexpr (dim > 1)
                v = v_half;
            else
                v = v_next += (v_half - v_next) * method.s_star_mult;

            if (x <= x_max) {
                soln.back().c_minus = minus_count;
                soln.back().c_plus = plus_count;
                soln.push_back({x,v_next,v_half,s,h,0,0});
            }
            minus_count = 0;
            plus_count = 0;
            if (std::abs(s) < eps_lower_bound) {
                h_half = h;
                h *= 2.0;
                ++plus_count;
            }
        }

        ++iteration;
    }

    for (size_t i = soln.size()-1; i > 0; --i) {
        soln[i].c_plus  = soln[i-1].c_plus;
        soln[i].c_minus = soln[i-1].c_minus;
    }

    return {soln, iteration};
}


////////////////////////////////////////////////////////////////////////////////////


// Реализации шагов методов

Point<2> rk4_system_2_next_point(double x, const Point<2> &v, const Diff_equation<Point<2>, 2> &F, double h) {
    Point<2> k1, k2, k3, k4;
    double h2 = 0.5 * h;

    k1[0] = F[0](x,v);
    k1[1] = F[1](x,v);
    k2[0] = F[0](x + h2, v + h2*k1);
    k2[1] = F[1](x + h2, v + h2*k1);
    k3[0] = F[0](x + h2, v + h2*k2);
    k3[1] = F[1](x + h2, v + h2*k2);
    k4[0] = F[0](x + h , v + h *k3);
    k4[1] = F[1](x + h , v + h *k3);

    return v + h/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
}

Method<Point<2>, 2> rk4_system2(rk4_system_2_next_point, 4); // Метод Рунге-Кутта 4 порядка для систем ОДУ 2 порядка
