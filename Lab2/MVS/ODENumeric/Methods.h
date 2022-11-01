#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include <array>
#include <limits>


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

// Расстояние между точками (нужно для получения оценки локальной погрешности)
template <size_t dim>
double dist(const Point<dim> &a, const Point<dim> &b) {
	double result = 0.0;
	for (size_t i = 0; i < dim; ++i)
		result += (b[i] - a[i]) * (b[i] - a[i]);
	return sqrt(result);
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
struct Ivp_solution {
	std::vector<double> x, h, s_star;
	std::vector<size_t> c_plus, c_minus;
	std::vector<value_t> v, v2;

	Ivp_solution(const std::vector<double> &x       = std::vector<double>(),
				 const std::vector<value_t> &v      = std::vector<value_t>(),
				 const std::vector<value_t> &v2     = std::vector<value_t>(),
				 const std::vector<double> &h       = std::vector<double>(),
				 const std::vector<double> &s_star  = std::vector<double>(),
				 const std::vector<size_t> &c_plus  = std::vector<size_t>(),
				 const std::vector<size_t> &c_minus = std::vector<size_t>()):
		x(x), v(v), v2(v2), h(h), s_star(s_star), c_plus(c_plus), c_minus(c_minus) {}
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
Ivp_solution<value_t, dim> ivp_no_step_adjust (const Method<value_t, dim> &method,
											 const Diff_equation<value_t, dim> &F,
											 const double x0, 
											 const value_t &u0, 
											 const double h, 
											 const double x_max,
											 const double x_eps)
{
	double x = x0;
	value_t v(u0);
	Ivp_solution<value_t, dim> soln;

	soln.x.push_back(x);
	soln.v.push_back(v);
	soln.h.push_back(h);
	while (x_max - x > x_eps) {
		v = method.next_point(x,v,F,h);
		x += h;

		soln.x.push_back(x);
		soln.v.push_back(v);
		soln.h.push_back(h);
	}

	return soln;
}


// Общий вид метода с регулировкой шага
template <typename value_t, size_t dim>
Ivp_solution<value_t, dim> ivp_step_adjust (const Method<value_t, dim> &method,
										  const Diff_equation<value_t, dim> &F,
										  const double x0, 
										  const value_t &u0, 
										  double h, 
										  const double x_max,
										  const double x_eps,
										  const double v_eps, 
										  int nmax = -1) 
{

	// Если максимальное число итераций не задано, вычисляем его
	if (nmax == -1)
		nmax = 100*(int)((x_max - x0)/h);
	int iteration = 0;
	int plus_count = 0;
	int minus_count = 0;

	double x = x0;
	value_t v = u0, v_next, v_half;
	double h_half = 0.5*h;
	double s;       // контрольная величина;
	double s_star;  // оценка локальной погрешности
	double eps_lower_bound = v_eps/pow(2.0, method.p + 1); // при |S| < eps/2^(p+1) увеличивается шаг

	Ivp_solution<value_t, dim> soln;
	soln.x.push_back(x);
	soln.v.push_back(v);
	soln.v2.push_back(v);
	soln.h.push_back(h);
	soln.s_star.push_back(0);
	soln.c_plus.push_back(0);
	soln.c_minus.push_back(0);

	while (iteration < nmax && x_max - x > x_eps) {
		v_half = method.next_point(x + h_half, method.next_point(x, v, F, h_half), F, h_half);
		v_next = method.next_point(x, v, F, h);
		s      = method.s_mult * dist(v_half, v_next);
		s_star = method.pow2 * s;

		if (std::abs(s) > v_eps || x_max - x + h < x_eps) {
			h = h_half;
			h_half *= 0.5;
			++plus_count;
		}
		else {
			x += h;
			if constexpr (dim > 1)
				v = v_half;
			else
				v = v_next + s_star;

			if (x <= x_max) {
				soln.x.push_back(x);
				soln.v.push_back(v);
				soln.v2.push_back(v_half);
				soln.h.push_back(h);
				soln.s_star.push_back(s_star);
				soln.c_plus.push_back(plus_count);
				soln.c_minus.push_back(minus_count);
			}
			if (std::abs(s) < eps_lower_bound) {
				h_half = h;
				h *= 2.0;
				++minus_count;
			}

			minus_count = 0;
			plus_count = 0;
		}

		++iteration;
	}

	return soln;
}


////////////////////////////////////////////////////////////////////////////////////


// Реализации шагов методов
double euler_next_point(double x, double v, const Diff_equation<double,1> &F, double h) {
	return v + h*F[0](x,v);
};

double rk4_next_point(double x, double v, const Diff_equation<double,1> &F, double h) {
	double k1, k2, k3, k4;
	double h2 = 0.5 * h;
	auto f = F[0];

	k1 = f(x,v);
	k2 = f(x + h2, v + h2*k1);
	k3 = f(x + h2, v + h2*k2);
	k4 = f(x + h,  v + h *k3);

	return v + h/6.0 * (k1 + 2.0*(k2 + k3) + k4);
};

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

Method<double, 1> rk4(rk4_next_point, 4); // Метод Рунге-Кутта 4 порядка для ОДУ
Method<Point<2>, 2> rk4_system2(rk4_system_2_next_point, 4); // Метод Рунге-Кутта 4 порядка для систему ОДУ 2 порядка


