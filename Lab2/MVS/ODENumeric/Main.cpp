#define _USE_MATH_DEFINES
#include "tridiagonal.h"
#include <iostream>
using namespace std;

//namespace plt = matplotlibcpp;

/*

auto f1 = [](double x, double u) {
	return 3*u;
};

auto f2 = [](double x, double u) {
	return 1/(x*x*x*x + 1)*u*u + u - u*u*u*sin(10.0*x);
};

Diff_equation<Point<2>, 2> lab_system{
	[](double x, const Point<2> &p) {
		return p[1];
	},
	[](double x, const Point<2> &p) {
		return -(p[1]*p[1] + sin(p[0]));
	}
};
Point<2> iv{-1.0, 0.2};

int main() {
	//auto sol2 = euler_ivp(0.0, 1.0, f, 0.01, 2.0);
	//auto sol = ivp_no_step_adjust<double, 1>(rk4, Diff_equation<double, 1>{f1}, 0.0, 2.0, 0.0001, 3.0, 1e-4);
	auto sol = ivp_step_adjust<double, 1>(rk4, Diff_equation<double, 1>{f2}, 0.4, 0.0, 0.01, 0.6, 1e-6, 1e-4, 10000);
	//auto sol = ivp_no_step_adjust<Point<2>, 2>(rk4_system2, lab_system, 0.0, iv, -0.001, 2.5, 1e-3);
	//auto sol = ivp_step_adjust<Point<2>, 2>(rk4_system2, lab_system, 0.0, iv, 0.001, 10, 1e-3, 1e-6, 1000);

	cout << "x           u\n";
	for (int i = 0; i < sol.x.size(); ++i) {
		cout.precision(8);
		cout << fixed << sol.x[i] << "  ";
		cout.precision(16);
		cout << sol.v[i] << '\n';
		//cout << sol.v[i][0] << "  " << sol.v[i][1] << '\n'; //<< exp(3.0*sol1.x[i]) << "  " << abs(sol1.v[i] - exp(3.0*sol1.x[i])) << '\n';
	}

	//plt::plot(sol1.x, sol1.v);
}
*/

const double C1 = 0.1082099375634999;
const double C2 = 1.8917900624365001;
const double C3=  0.1427998568251809;
const double C4 = 2.1641463834260162;
const double C5 = M_PI/(2.0*sqrt(2.0));
const double C6 = 8*sqrt(2)/(M_PI * M_PI);
const double INV_3 = 1/3.0;

inline double u1(double x) {
	return -1.0 + C1*exp(x) + C2*exp(-x);
}

inline double u2(double x) {
	return -C6 + C3*exp(C5*x) + C4*exp(-C5*x);
}

inline double k1(double a, double b) {
	return 1/sqrt(2)*(log(tan(b/2.0 + 0.1)/tan(a/2.0 + 0.1)));
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
	
	const size_t n = 1000;
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
	v =  solve(n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	v2 = solve(2*n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	*/
	v = solve(n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1t, k2t, q1t, q2t, f1t, f2t);
	cout.precision(16);
	for (size_t i = 0; i <= n; ++i)
		cout << fixed << x[i] << "  " << u[i] << "  " << v[i] << "  " << abs(u[i] - v[i]) << '\n';
		//cout << fixed << x[i] << "  " << v[i] << "  " << v2[2*i] << "  " << abs(v[i] - v2[2*i]) << '\n';

	return 0;
}

