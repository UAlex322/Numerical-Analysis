#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include <array>
#include <limits>
#include <iostream>


using func = std::function<double(double, double)>;
std::vector<double> solve(size_t n, double a, double b, double m1, double m2, double xi, func k1, func k2, func q1, func q2, func f1, func f2) {
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