#ifndef WORKER_H
#define WORKER_H

#include <QObject>

const double a = 0.0, b = 2.0, c = 0.0, d = 1.0;

struct solution_t {
    size_t n, m;
    size_t k;
    std::vector<size_t> row_bound;  // границы строк, входящих в сетку
    std::vector<double> vVec;       // вектор численного решения
    double r_max = 0.0;             // максимум-норма невязки на выходе
    double r0_max = 0.0;            // максимум-норма невязки в начале
    double r_2 = 0.0;               // евклидова норма невязки на выходе
    double eps;                     // точность метода на выходе
    size_t iteration_num;           // число итераций на выходе
    double eps_max;
    size_t iteration_num_max;
};

Q_DECLARE_METATYPE(solution_t);

class Worker: public QObject {
    Q_OBJECT

public:
    explicit Worker(QObject* parent = 0);
    solution_t solve_sle_chebyshev_custom_region(
                double a, double b, double c, double d,
                size_t n, size_t m,
                const std::function<double(double, double)> &mu,
                const std::function<double(double, double)> &F,
                size_t iterations, double epsilon,
                size_t K, double Mmax, double Mmin);
public slots:
    void launchMethod(size_t n, size_t m, size_t K, double eps, size_t nmax);
signals:
    void updateProgress(const QString&);
    void sendResult(const solution_t&);
    void finished();
};

double u(double x, double y);

double f(double x, double y);

#endif // WORKER_H
