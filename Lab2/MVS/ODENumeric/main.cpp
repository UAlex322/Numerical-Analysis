#include "tridiagonal.h"
#include <iostream>
using namespace std;

int main() {
	
	const size_t n = 20000;
	vector<double> x(n+1), u(n+1), v, v2;
	double xc = 0.0;
	double h = 1.0/(double)n;

	u = get_true_test_solution(n);
	for (size_t i = 0; i <= n; ++i, xc += h)
		x[i] = xc;
	
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

	/*
	v =  solve(n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	v2 = solve(2*n, 0.0, 1.0, 1.0, 0.0, M_PI_4, k1, k2, q1, q2, f1, f2);
	double max_err = 0.0;
	for (size_t i = 0; i <= n; ++i)
		if (abs(v[i] - v2[i*2]) > max_err)
			max_err = abs(v[i] - v2[i*2]);
	cout.precision(16);
	cout << "Main problem\n";
	cout << "N = " << n << "\nMaximal error: " << max_err;
	*/
	//for (size_t i = 0; i <= n; ++i)
		//cout << fixed << x[i] << "  " << u[i] << "  " << v[i] << "  " << abs(u[i] - v[i]) << '\n';
		//cout << fixed << x[i] << "  " << v[i] << "  " << v2[2*i] << "  " << abs(v[i] - v2[2*i]) << '\n';

	return 0;
}

