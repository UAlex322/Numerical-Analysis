#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include <array>
#include <limits>
#include <iostream>


void solveMatrix(int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
	{
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
	}
}

using func = std::function<double(double, double)>;
std::vector<double> solve(size_t n, double a, double b, double m1, double m2, double xi, func k1, func k2, func q1, func q2, func f1, func f2) {
	double h = (b-a)/n;                                   // шаг численного метогда
	size_t discont_index_a  = (size_t)(xi/h) + 1;         // индекс, на котором происходит разрыв функций при движении по точкам a + h*i
	size_t discont_index_dp = (size_t)(xi/h - 0.5) + 1;   // индекс, на котором происходит разрыв функций при движении по точкам a + h*(i + 1/2)
	std::vector<double> ai(n+1), di(n), pi(n);            // массивы коэффициентов дл¤ построени¤ системы

	// ¬џ„»—Ћ≈Ќ»≈  ќЁ‘‘»÷»≈Ќ“ќ¬ —Ћј”
	
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
	


	/*
	for (size_t i = 1; i < discont_index_a; ++i)
		ai[i] = h/k1(a + (i-1)*h, a + i*h);
	for (size_t i = discont_index_a + 1; i <= n; ++i)
		ai[i] = h/k2(a + (i-1)*h, a + i*h);
	for (size_t i = 1; i < discont_index_dp; ++i) {
		di[i] = q1(a + h/2 + (i-1)*h, a + h/2 + i*h)/h;
		pi[i] = f1(a + h/2 + (i-1)*h, a + h/2 + i*h)/h;
	}
	for (size_t i = discont_index_dp + 1; i < n; ++i) {
		di[i] = q2(a + h/2 + (i-1)*h, a + h/2 + i*h)/h;
		pi[i] = f2(a + h/2 + (i-1)*h, a + h/2 + i*h)/h;
	}
	ai[discont_index_a] = h/(k1(a + (discont_index_a-1)*h, xi) + k2(xi, a + discont_index_a*h));
	di[discont_index_dp] = (q1(a + h/2 + (discont_index_dp-1)*h, xi) + q2(xi, a + h/2 + discont_index_dp*h))/h;
	pi[discont_index_dp] = (f1(a + h/2 + (discont_index_dp-1)*h, xi) + f2(xi, a + h/2 + discont_index_dp*h))/h;
	*/


	double hi = 1/h;                      // обратный шаг (множитель в выражени¤х ниже)
	std::vector<double> A(n), B(n), C(n); // коэффициенты на диагонал¤х матрицы системы
	for (size_t i = 1; i < n; ++i) {
		A[i] = ai[i] * hi;                             // ai[i]/(h^2)
		B[i] = ai[i+1] * hi;                           // ai[i+1]/(h^2)
		C[i] = (ai[i] + ai[i+1]) * hi + di[i] * h; // (ai[i] + ai[i+1])/(h^2) + di[i]
		pi[i] *= h;
	}

	// ћ≈“ќƒ ѕ–ќ√ќЌ »

	std::vector<double> alpha(n+1), beta(n+1); // коэффициенты метода прогонки
	std::vector<double> y(n+1);                // вектор решени¤

	
	alpha[1] = 0, beta[1] = m1;                // начальные значени¤ коэффициентов метода прогонки
	for (size_t i = 1; i < n; ++i) {
		double divide = 1/(C[i] - alpha[i]*A[i]); // знаменатель в коэффициентах метода прогонки
		alpha[i+1] = B[i] * divide;
		beta[i+1]  = (pi[i] + beta[i]*A[i]) * divide;
	}
	//y[0] = m1; // на основании начальных условий
	y[n] = m2; // на основании начальных условий
	for (size_t i = n; i >= 1; --i)
		y[i-1] = alpha[i]*y[i] + beta[i];


	/*
	std::vector<double> A_(n+1), B_(n+1), C_(n+1), F_(n+1), y_(n+1);
	for (int i = 1; i < n; ++i)
		A_[i] = A[i];
	A_[n] = 0.0;
	for (int i = 1; i < n; ++i)
		B_[i] = B[i];
	B_[0] = 0.0;
	for (int i = 1; i < n; ++i)
		C_[i] = -C[i];
	C_[0] = C_[n] = 1.0;
	for (int i = 1; i < n; ++i)
		F_[i] = -pi[i];
	F_[0] = m1;
	F_[n] = m2;
	solveMatrix(n+1, A_.data(), C_.data(), B_.data(), F_.data(), y_.data());
	


	std::cout << "COMPARISON\n";
	for (int i = 0; i <= n; ++i) {
		std::cout << y[i] << "  " << y_[i] << '\n';
	}
	*/
	return y;
}

std::vector<double> tridiagonal_solve(size_t n,
									  const std::vector<double> &a,
									  const std::vector<double> &b,
									  const std::vector<double> &c,
									  const std::vector<double> &p,
									  double k1,
									  double k2,
									  double m1,
									  double m2) {
	std::vector<double> y(n+1), alpha(n+1), beta(n+1);

	alpha[1] = k1, beta[1] = m1;
	for (size_t i = 1; i < n; ++i) {
		double divide = 1/(c[i] - alpha[i]*a[i]);
		alpha[i+1] =  b[i] * divide;
		beta[i+1]  = (p[i] + beta[i]*a[i]) * divide;
	}
	y[n] = (m2 + k2*beta[n])/(1 - k2*alpha[n]);
	for (size_t i = n; i >= 1; --i)
		y[i-1] = alpha[i]*y[i] + beta[i];

	return y;
}