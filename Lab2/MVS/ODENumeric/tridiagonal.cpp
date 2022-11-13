#define _USE_MATH_DEFINES
#include <functional>
#include <cmath>
#include <array>
#include <limits>
#include "tridiagonal.h"


using func = std::function<double(double, double)>;
std::vector<double> solve_heat_equation(size_t n, double a, double b, double m1, double m2, double xi, func k1, func k2, func q1, func q2, func f1, func f2) {
	double h = (b-a)/n;                                   // шаг численного метогда
	size_t discont_index_a  = (size_t)(xi/h) + 1;         // индекс, на котором происходит разрыв функций при движении по точкам a + h*i
	size_t discont_index_dp = (size_t)(xi/h - 0.5) + 1;   // индекс, на котором происходит разрыв функций при движении по точкам a + h*(i + 1/2)
	std::vector<double> ai(n+1), di(n), pi(n);            // массивы коэффициентов дл¤ построени¤ системы

	double x = a; // текуща¤ точка вида 'x + h*i'
	for (size_t i = 1; i < discont_index_a; ++i, x += h) {
		ai[i] = h/k1(x, x+h); // вычисление коэффициента по первой функции
	}
	ai[discont_index_a] = h/(k1(x, xi) + k2(xi, x+h)); // отдельное вычисление коэффициента на отрезке с точкой разрыва
	x += h;
	for (size_t i = discont_index_a + 1; i <= n; ++i, x += h) {
		ai[i] = h/k2(x, x+h); // вычисление коэффициента по второй функции
	}

	x = a + h/2; // текуща¤ точка вида 'x + h*(i + 1/2)'
	for (size_t i = 1; i < discont_index_dp; ++i, x += h) {
		di[i] = q1(x, x+h)/h; // вычисление коэффициента по первой функции
		pi[i] = f1(x, x+h)/h; // вычисление коэффициента по первой функции
	}
	di[discont_index_dp] = (q1(x, xi) + q2(xi, x+h))/h; // отдельное вычисление коэффициента на отрезке с точкой разрыва
	pi[discont_index_dp] = (f1(x, xi) + f2(xi, x+h))/h; // отдельное вычисление коэффициента на отрезке с точкой разрыва
	x += h;
	for (size_t i = discont_index_dp + 1; i < n; ++i, x += h) {
		di[i] = q2(x, x+h)/h; // вычисление коэффициента по второй функции
		pi[i] = f2(x, x+h)/h; // вычисление коэффициента по второй функции
	}


	double hi = 1/h;                      // обратный шаг (множитель в выражени¤х ниже)
	std::vector<double> A(n), B(n), C(n); // коэффициенты на диагонал¤х матрицы системы
	for (size_t i = 1; i < n; ++i) {
		A[i] = ai[i] * hi;                         // ai[i]/(h^2) * h
		B[i] = ai[i+1] * hi;                       // ai[i+1]/(h^2) * h
		C[i] = (ai[i] + ai[i+1]) * hi + di[i] * h; // ((ai[i] + ai[i+1])/(h^2) + di[i]) * h
		pi[i] *= h;
	}


	std::vector<double> alpha(n+1), beta(n+1); // коэффициенты метода прогонки
	std::vector<double> y(n+1);                // вектор решени¤
	alpha[1] = 0, beta[1] = m1;                // начальные значени¤ коэффициентов метода прогонки
	for (size_t i = 1; i < n; ++i) {
		double divide = 1/(C[i] - alpha[i]*A[i]); // знаменатель в коэффициентах метода прогонки
		alpha[i+1] = B[i] * divide;
		beta[i+1]  = (pi[i] + beta[i]*A[i]) * divide;
	}
	y[n] = m2; // на основании начальных условий
	for (size_t i = n; i >= 1; --i)
		y[i-1] = alpha[i]*y[i] + beta[i];

	return y;
}


// Константы, используемые в функциях

const double C1 = -0.3393176035227834;
const double C2 = 0.3393176035227834;
const double C3=  -0.4920418012319509;
const double C4 = 1.0560782612867800;
const double C5 = M_PI/(2.0*sqrt(2.0));
const double C6 = 8*sqrt(2)/(M_PI * M_PI);
const double INV_3 = 1/3.0;

// Функции, соответствующие участкам непрерывности истинного решения

inline double u1(double x) {
	return 1.0 + C1*exp(x) + C2*exp(-x);
}

inline double u2(double x) {
	return C6 + C3*exp(C5*x) + C4*exp(-C5*x);
}

std::vector<double> get_true_test_solution(size_t n) {
	std::vector<double> u(n+1);
	double xc = 0.0;
	double h = 1.0/(double)n;
	const size_t discont_pos = (size_t)(M_PI_4 / h);

	xc = 0.0;
	for (size_t i = 0; i < discont_pos; ++i, xc += h)
		u[i] = u1(xc);
	for (size_t i = discont_pos; i <= n; ++i, xc += h)
		u[i] = u2(xc);

	return u;
}

// Функции, участвующие в решении основной задачи

inline double k1(double a, double b) {
	return M_SQRT1_2*(log(tan(b/2.0 + 0.1)/tan(a/2.0 + 0.1)));
}

inline double k2(double a, double b) {
	return tan(b) - tan(a);
}

inline double q1(double a, double b) {
	return b - a;
}

inline double q2(double a, double b) {
	return INV_3 * (b*b*b - a*a*a);
}

inline double f1(double a, double b) {
	return 0.5 * (cos(2.0*a) - cos(2.0*b));
}

inline double f2(double a, double b) {
	return sin(b) - sin(a);
}

std::vector<double> solve_main(size_t n) {
	return solve_heat_equation(n, 0.0, 1.0, 1.0, 0.0, M_PI/4.0, k1, k2, q1, q2, f1, f2);
}

// Функции, участвующие в решении тестовой задачи

inline double k1t(double a, double b) {
	return b-a;
}

inline double k2t(double a, double b) {
	return 2.0*(b-a);
}

inline double q1t(double a, double b) {
	return b-a;
}

inline double q2t(double a, double b) {
	return M_PI*M_PI/16.0*(b-a);
}

inline double f1t(double a, double b) {
	return b-a;
}

inline double f2t(double a, double b) {
	return (b-a)/sqrt(2);
}

std::vector<double> solve_test(size_t n) {
	return solve_heat_equation(n, 0.0, 1.0, 1.0, 0.0, M_PI/4.0, k1t, k2t, q1t, q2t, f1t, f2t);
}