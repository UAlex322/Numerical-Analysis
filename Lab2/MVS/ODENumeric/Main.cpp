#define _USE_MATH_DEFINES
#include "tridiagonal.h"
#include <iostream>
using namespace std;

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


int main() {
	
	const size_t n = 20000;
	vector<double> x(n+1), u(n+1), v, v2;
	double xc = 0.0;
	double h = 1.0/(double)n;
	const size_t discont_pos = (size_t)(M_PI_4 / h);

	for (size_t i = 0; i <= n; ++i, xc += h)
		x[i] = xc;
	xc = 0.0;
	for (size_t i = 0; i < discont_pos; ++i, xc += h)
		u[i] = u1(xc);
	for (size_t i = discont_pos; i <= n; ++i, xc += h)
		u[i] = u2(xc);
	
	/*
	v = solve(n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1t, k2t, q1t, q2t, f1t, f2t);
	double max_err = 0.0;
	for (size_t i = 0; i <= n; ++i)
		if (abs(u[i] - v[i]) > max_err)
			max_err = abs(u[i] - v[i]);
	cout.precision(16);
	cout << "Test problem\n";
	cout << "N = " << n << "\nMaximal error: " << max_err;
	*/

	
	v =  solve(n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	v2 = solve(2*n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	double max_err = 0.0;
	for (size_t i = 0; i <= n; ++i)
		if (abs(v[i] - v2[i*2]) > max_err)
			max_err = abs(v[i] - v2[i*2]);
	cout.precision(16);
	cout << "Main problem\n";
	cout << "N = " << n << "\nMaximal error: " << max_err;
	//for (size_t i = 0; i <= n; ++i)
		//cout << fixed << x[i] << "  " << u[i] << "  " << v[i] << "  " << abs(u[i] - v[i]) << '\n';
		//cout << fixed << x[i] << "  " << v[i] << "  " << v2[2*i] << "  " << abs(v[i] - v2[2*i]) << '\n';

	return 0;
}

