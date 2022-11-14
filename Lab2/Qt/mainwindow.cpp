#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QDebug>
#include <tridiagonal.h>
#include <qcustomplot.h>
#include <limits>
#include <QtMath>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->testPlot->setInteraction(QCP::iRangeZoom,true);
    ui->testPlot->setInteraction(QCP::iRangeDrag,true);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// ТЕСТОВОЕ ЗАДАНИЕ
void MainWindow::on_pushButton_pressed()
{
    QString str_n = ui->lineEdit->text();

    size_t n;

    if (!checkInput(str_n, n))
    {
        return;
    }

    auto stdv = solve_test(n);
    QVector<double> v(stdv.begin(), stdv.end());
    auto stdu = solve_test(n);
    QVector<double> u(stdu.begin(), stdu.end());

    QVector<double> x(n+1), xi(n+1);
    double h = 1.0/(double)n;
    double xc = 0.0;
    for (size_t i = 0; i <= n; ++i, xc += h)
        x[i] = xc;
    for (size_t i = 0; i <= n; ++i)
        xi[i] = u[i] - v[i];

    auto max_err_iter = std::max_element(xi.begin(), xi.end());
    double max_err = *max_err_iter;
    double max_err_pos = (max_err_iter - xi.begin())*h;


    // График 1
    QCustomPlot* customPlot1 = ui->testPlot;
    customPlot1->clearGraphs();

    customPlot1->xAxis->setRange(0, 1);
    customPlot1->yAxis->setRange(0, 1);

    customPlot1->xAxis->setLabel("x");
    customPlot1->yAxis->setLabel("u");

    customPlot1->addGraph();
    customPlot1->graph(0)->addData(x, u);
    customPlot1->graph(0)->setPen(QPen(Qt::blue));
    customPlot1->graph(0)->setName("Аналитическое решение u(x)");

    customPlot1->addGraph();
    customPlot1->graph(1)->addData(x, v);
    customPlot1->graph(1)->setPen(QPen(Qt::green));
    customPlot1->graph(1)->setName("Численное решение v(x)");

    customPlot1->addGraph();
    customPlot1->graph(2)->addData(x, xi);
    customPlot1->graph(2)->setPen(QPen(Qt::red));
    customPlot1->graph(2)->setName("Разность аналитического и численного решения |u(x)-v(x)|");

    customPlot1->replot();
}

// ОСНОВНОЕ ЗАДАНИЕ
void MainWindow::on_main1SolveButton_clicked()
{
    QString str_n = ui->lineEdit->text();

    size_t n;
    if (!checkInput(str_n, n))
    {
        return;
    }

    auto stdv = solve_test(n);
    QVector<double> v(stdv.begin(), stdv.end());
    auto stdv2 = solve_test(n);
    QVector<double> v2(stdv2.begin(), stdv2.end());

    QVector<double> x(n+1), xi(n+1);
    double h = 1.0/(double)n;
    double xc = 0.0;
    for (size_t i = 0; i <= n; ++i, xc += h)
        x[i] = xc;
    for (size_t i = 0; i <= n; ++i)
        xi[i] = v[i] - v2[i];

    auto max_err_iter = std::max_element(xi.begin(), xi.end());
    double max_err = *max_err_iter;
    double max_err_pos = (max_err_iter - xi.begin())*h;

    // График 1
    QCustomPlot* customPlot1 = ui->testPlot;
    customPlot1->clearGraphs();

    customPlot1->xAxis->setRange(0, 1);
    customPlot1->yAxis->setRange(-0.5, 1);

    customPlot1->xAxis->setLabel("x");
    customPlot1->yAxis->setLabel("v");

    customPlot1->addGraph();
    customPlot1->graph(0)->addData(x, v);
    customPlot1->graph(0)->setPen(QPen(Qt::blue));
    customPlot1->graph(0)->setName("Аналитическое решение u(x)");

    customPlot1->addGraph();
    customPlot1->graph(1)->addData(x, v2);
    customPlot1->graph(1)->setPen(QPen(Qt::green));
    customPlot1->graph(1)->setName("Численное решение v(x)");

    customPlot1->addGraph();
    customPlot1->graph(2)->addData(x, xi);
    customPlot1->graph(2)->setPen(QPen(Qt::red));
    customPlot1->graph(2)->setName("Разность аналитического и численного решения |v(x)-v2(x)|");

    customPlot1->replot();
}

bool MainWindow::checkInput(const QString& str_n, size_t& n)
{
    bool result = false;

    if (str_n.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'n' должно быть непустым");
    else
        result = true;

    if (result == false)
        return false;
    n = str_n.toInt(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'n' имеет некорректное значение");
        return false;
    }

    result = false;
    if (n <= 0)
        QMessageBox::critical(this, "Ошибка!", "Число разбиений 'n' должен быть положительным");
    else
        result = true;

    return result;
}
