#include "dirichlet.h"
#include <clocale>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
using namespace std;

template <typename ValueT, typename HandlerT>
void checkInput(ValueT &value, const string &outText, const string &errorText, HandlerT &handler) {
    bool success = false;
    string inputString;
    do {
        try {
            cout << outText;
            cin >> inputString;
            value = handler(inputString);

            success = true;
        }
        catch(...) {
            cout << errorText << '\n';
        }
    } while (!success);
}

double doubleHandler(string &inputString) {
    double val = stof(inputString);
    if (val <= 0.0)
        throw;
    return val;
}

int64_t intHandler(string &inputString) {
    int64_t val = stoll(inputString);
    if (val <= 0)
        throw;
    return val;
}

void print_solution(const vector<double> &u, size_t n, size_t m, size_t precision) {
    double x;
    double y;
    double h = 1.0/n;
    double k = 1.0/m;
    size_t width = precision + 4;

    cout.precision(precision);

    y = 1.0;
    for (int64_t j = m; j >= 0; --j, y -= k) {
        x = 0.0;
        cout << setw(width) << fixed << y;
        for (size_t i = 0; i <= n; ++i, x += h)
            cout << setw(width) << fixed << u[(n+1)*j + i];
        cout << '\n';
    }
    cout << setw(width) << "y/x";
    for (size_t i = 0; i <= n; ++i)
        cout << setw(width) << fixed << h*i;
    cout << "\n\n";
}

const string introMessage = "\
\"��������� ������� �������� ������ ������� ��� ��������� �������� �� ������������ �������\".\n\
������ 382003_3, ������� 5 (�����, ��������, �������, �������).\n\n\
�������: ������������� [0; 1] x [0; 1].\n\n\
������ �������:\n\
    u''xx + u''yy = 6x + 2, ���� (x,y) ������ ��������������;\n\
    u(x,y) = x^3 + y^2 + 3, ���� (x,y) ����� �� �������.\n\
������ �������: u(x,y) = x^3 + y^2 + 3.\n\n\
��������� ������� ������: � ������� ���������� �����.\n\
�����: ����������� �� 'x' (� ����� h = 1/n) � 'y' (� ����� k = 1/m),\n\
    n > 0, m > 0 - �������� ������������� ��������� �����.\n\
����� ������� ���� ���������� �����: ������������ (�������).\n\
�������� ���������:\n\
a) ���������� ��������� ������������� ����� �������� (N_max > 0);\n\
�) ���������� �������� ������������� �������� (eps_max > 0):\n\
   ||v^(s) - v^(s-1)|| < eps_max (s - ����� �������� ���� ������, (s-1) - ����� ����������� ����).\n\
���������: ������ ������� ������ ������� ��������� � ������ �������� ���������� �����,\n\
������� ����������� ����� ||Z|| ��������� � ������������ ������������� ������.\n\n\
";

int main() {
    setlocale(LC_ALL, "Russian");

    auto U = [](double x, double y) -> double {
        return x*x*x + y*y + 3.0;
    };
    auto F = [](double x, double y) -> double {
        return -6.0*x - 2.0;
    };
    double a = 0.0;
    double b = 1.0;
    double c = 0.0;
    double d = 1.0;
    int64_t n;
    int64_t m;
    int64_t N_max;
    int64_t N = 0;
    double eps = 0.0;
    double eps_max;
    int64_t prec;

    cout << introMessage;
    
    while (true) {
        checkInput(n, "������� �������� ����� (n): ", "������������ �������� n, ��������� ����", intHandler);
        checkInput(m, "������� �������� ����� (m): ", "������������ �������� m, ��������� ����", intHandler);
        checkInput(N_max, "������� ������������ ����� ����� (N_max): ", "������������ �������� N_max, ��������� ����", intHandler);
        checkInput(eps_max, "������� �������� (eps_max): ", "������������ �������� eps_max, ��������� ����", doubleHandler);
        checkInput(prec, "������� ����� ������ ����� ������� � ���������� �������: ", "������������ �������� ����� ������, ��������� ����", intHandler);

        auto v = solve_dirichlet_optimized(a, b, c, d, n, m, U, U, U, U, F, N_max, N, eps_max, eps);
        cout << "\n���������:\n";
        cout << "��������� N = " << N << " ���(��).\n";
        cout << "�������� �� ������ eps_N = " << eps << '\n';
        cout << "��������� ����� ������� �� ������ ||R_N|| = " << discrepancy(F, v, n, m) << '\n';
        cout << "��������-����� ����������� �����  ||Z_N|| = " << solution_error(U, v, n, m) << '\n';
        cout << '\n';

        vector<double> u((n+1)*(m+1));
        for (size_t j = 0; j <= m; ++j)
            for (size_t i = 0; i <= n; ++i)
                u[(n+1)*j + i] = U((double)i/n, (double)j/m);
        cout << "��������� �������:\n";
        print_solution(v, n, m, prec);
        cout << "\n������ �������:\n";
        print_solution(u, n, m, prec);
    }

    return 0;
}

