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
\"Численное решение тестовой задачи Дирихле для уравнения Пуассона на прямугольной области\".\n\
Группа 382003_3, команда 5 (Ларин, Медведев, Устинов, Шипилов).\n\n\
Область: прямоугольник [0; 1] x [0; 1].\n\n\
Задача Дирихле:\n\
    u''xx + u''yy = 6x + 2, если (x,y) внутри прямоугольника;\n\
    u(x,y) = x^3 + y^2 + 3, если (x,y) лежит на границе.\n\
Точное решение: u(x,y) = x^3 + y^2 + 3.\n\n\
Численное решение задачи: с помощью разностной схемы.\n\
Сетка: равномерная по 'x' (с шагом h = 1/n) и 'y' (с шагом k = 1/m),\n\
    n > 0, m > 0 - заданные пользователем параметры сетки.\n\
Метод решения СЛАУ разностной схемы: итерационный (Зейделя).\n\
Критерии остановки:\n\
a) выполнение заданного пользователем числа итераций (N_max > 0);\n\
б) достижение заданной пользователем точности (eps_max > 0):\n\
   ||v^(s) - v^(s-1)|| < eps_max (s - номер текущего шага метода, (s-1) - номер предыдущего шага).\n\
ЗАМЕЧАНИЕ: точное решение задачи Дирихле совпадает с точным решением разностной схемы,\n\
поэтому погрешность схемы ||Z|| совпадает с погрешностью итерационного метода.\n\n\
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
        checkInput(n, "Введите параметр сетки (n): ", "Некорректное значение n, повторите ввод", intHandler);
        checkInput(m, "Введите параметр сетки (m): ", "Некорректное значение m, повторите ввод", intHandler);
        checkInput(N_max, "Введите максимальное число шагов (N_max): ", "Некорректное значение N_max, повторите ввод", intHandler);
        checkInput(eps_max, "Введите точность (eps_max): ", "Некорректное значение eps_max, повторите ввод", doubleHandler);
        checkInput(prec, "Введите число знаков после запятой в выведенном решении: ", "Некорректное значение числа знаков, повторите ввод", intHandler);

        auto v = solve_dirichlet_optimized(a, b, c, d, n, m, U, U, U, U, F, N_max, N, eps_max, eps);
        cout << "\nРЕЗУЛЬТАТ:\n";
        cout << "Выполнено N = " << N << " шаг(ов).\n";
        cout << "Точность на выходе eps_N = " << eps << '\n';
        cout << "Евклидова норма невязки на выходе ||R_N|| = " << discrepancy(F, v, n, m) << '\n';
        cout << "Максимум-норма погрешности схемы  ||Z_N|| = " << solution_error(U, v, n, m) << '\n';
        cout << '\n';

        vector<double> u((n+1)*(m+1));
        for (size_t j = 0; j <= m; ++j)
            for (size_t i = 0; i <= n; ++i)
                u[(n+1)*j + i] = U((double)i/n, (double)j/m);
        cout << "Численное решение:\n";
        print_solution(v, n, m, prec);
        cout << "\nТочное решение:\n";
        print_solution(u, n, m, prec);
    }

    return 0;
}

