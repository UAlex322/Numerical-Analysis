#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QDebug>
#include <tridiagonal.h>
#include <qcustomplot.h>
#include <limits>
#include <QtMath>

void MainWindow::InfoWidget(QWidget *parent, QLabel *label)
{
    parent->setWindowTitle("Справка о решении");
    QHBoxLayout *layout = new QHBoxLayout(parent);
    layout->addWidget(label);
    parent->setLayout(layout);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    window = new QWidget;
    windowLabel = new QLabel(window);
    InfoWidget(window, windowLabel);
    window->activateWindow();
    setAttribute(Qt::WA_DeleteOnClose);

    ui->setupUi(this);

    ui->tPlot->setInteraction(QCP::iRangeZoom,true);
    ui->tPlot->setInteraction(QCP::iRangeDrag,true);

    ui->tTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->tTable->setStyleSheet("border: 1px solid grey");
    ui->tTable->setModel(new QStandardItemModel(0,4,this));
    ui->tTable->model()->setHeaderData(0, Qt::Horizontal, "x_i");
    ui->tTable->model()->setHeaderData(1, Qt::Horizontal, "u_i");
    ui->tTable->model()->setHeaderData(2, Qt::Horizontal, "v_i");
    ui->tTable->model()->setHeaderData(3, Qt::Horizontal, "|u_i - v_i|");
    ui->tTable->setColumnWidth(0,75);
    for (int i = 1; i < 4; ++i)
        ui->tTable->setColumnWidth(i,90);


    ui->mPlot->setInteraction(QCP::iRangeZoom,true);
    ui->mPlot->setInteraction(QCP::iRangeDrag,true);

    ui->mTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->mTable->setStyleSheet("border: 1px solid grey");
    ui->mTable->setModel(new QStandardItemModel(0,4,this));
    ui->mTable->model()->setHeaderData(0, Qt::Horizontal, "x_i");
    ui->mTable->model()->setHeaderData(1, Qt::Horizontal, "v_i");
    ui->mTable->model()->setHeaderData(2, Qt::Horizontal, "v2_i");
    ui->mTable->model()->setHeaderData(3, Qt::Horizontal, "|v_i - v2_i|");
    ui->mTable->setColumnWidth(0,75);
    for (int i = 1; i < 4; ++i)
        ui->mTable->setColumnWidth(i,90);

    ui-> tPicN->setGeometry(145, 210, 16, 16);
    ui->tLineN->setGeometry(165, 210, 81, 16);
    ui-> mPicN->setGeometry(145, 210, 16, 16);
    ui->mLineN->setGeometry(165, 210, 81, 16);
    ui->tButtonSolnInfo->setGeometry(1100, 235, 141, 24);
    ui->mButtonSolnInfo->setGeometry(1100, 235, 141, 24);

    ui->tabWidget->setCurrentIndex(0);
}

MainWindow::~MainWindow()
{
    window->close();
    delete ui;
}

bool MainWindow::checkInput(const QString& str_n, int& n)
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

// ТЕСТОВОЕ ЗАДАНИЕ
void MainWindow::on_tButtonSolve_pressed()
{
    QString str_n = ui->tLineN->text();

    int n;

    if (!checkInput(str_n, n))
    {
        return;
    }

    auto stdv = solve_test(n);
    QVector<double> v(stdv.begin(), stdv.end());
    auto stdu = get_true_test_solution(n);
    QVector<double> u(stdu.begin(), stdu.end());

    QVector<double> x(n+1), z(n+1);
    double h = 1.0/(double)n;
    double xc = 0.0;
    for (int i = 0; i <= n; ++i, xc += h)
        x[i] = xc;
    for (int i = 0; i <= n; ++i)
        z[i] = abs(u[i] - v[i]);

    auto model = (QStandardItemModel*)ui->tTable->model();
    model->setRowCount(n);
    for (int i = 0; i <= n; ++i)
    {
        model->setData(model->index(i,0), x[i]);
        model->setData(model->index(i,1), u[i]);
        model->setData(model->index(i,2), v[i]);
        model->setData(model->index(i,3), z[i]);
    }

    auto max_err_iter = std::max_element(z.begin(), z.end());
    double max_err = *max_err_iter;
    double max_err_pos = (max_err_iter - z.begin())*h;

    // График 1
    QCustomPlot* plot = ui->tPlot;
    plot->clearGraphs();

    plot->xAxis->setRange(-0.2, 1);
    plot->yAxis->setRange(-0.2, 1);

    plot->xAxis->setLabel("x");
    plot->yAxis->setLabel("u");
    plot->legend->setVisible(true);

    plot->addGraph();
    plot->addGraph();
    plot->addGraph();

    plot->graph(1)->addData(x, v);
    plot->graph(1)->setPen(QPen(Qt::green));
    plot->graph(1)->setName("Численное решение v(x)");

    plot->graph(2)->addData(x, z);
    plot->graph(2)->setPen(QPen(Qt::red));
    plot->graph(2)->setName("Модуль разности аналитического и численного решений |u(x)-v(x)|");

    n = std::min(50*n, 1000000);
    stdu = get_true_test_solution(n);
    u.resize(n+1);
    std::copy(stdu.begin(), stdu.end(), u.begin());
    x.resize(n+1);
    xc = 0.0;
    h = 1/(double)n;
    for (int i = 0; i <= n; ++i, xc += h)
        x[i] = xc;

    plot->graph(0)->addData(x, u);
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setName("Аналитическое решение u(x)");

    plot->replot();

    if (ui->tButtonSolnInfo->isVisible())
        ui->tButtonSolnInfo->close();
    windowLabel->setText("Для решения задачи использована равномерная сетка с числом разбиейний n ="+ QString::number(n) +".\n"
                                 "Задача должна быть решена с погрешностью не более Ɛ = 0.5*10^(-6).\n"
                                 "Задача решена с погрешностью Ɛ₁ = " + QString::number(max_err) + ".\n"
                                 "Максимальное отклонение аналитического и численного решения наблюдается в точке x = "+ QString::number(max_err_pos)+".");
    ui->tButtonSolnInfo->click();
}


// ОСНОВНОЕ ЗАДАНИЕ
void MainWindow::on_mButtonSolve_pressed()
{
    QString str_n = ui->mLineN->text();

    int n;

    if (!checkInput(str_n, n))
    {
        return;
    }

    auto stdv = solve_main(n);
    QVector<double> v(stdv.begin(), stdv.end());
    auto stdv2 = solve_main(2*n);
    QVector<double> v2(n+1);
    for (int i = 0; i <= n; ++i)
        v2[i] = stdv2[2*i];

    QVector<double> x(n+1), z(n+1);
    double h = 1.0/(double)n;
    double xc = 0.0;
    for (int i = 0; i <= n; ++i, xc += h)
        x[i] = xc;
    for (int i = 0; i <= n; ++i)
        z[i] = abs(v[i] - v2[i]);

    auto model = (QStandardItemModel*)ui->mTable->model();
    model->setRowCount(n);
    for (int i = 0; i <= n; ++i)
    {
        model->setData(model->index(i,0), x[i]);
        model->setData(model->index(i,1), v[i]);
        model->setData(model->index(i,2),v2[i]);
        model->setData(model->index(i,3), z[i]);
    }

    auto max_err_iter = std::max_element(z.begin(), z.end());
    double max_err = *max_err_iter;
    double max_err_pos = (max_err_iter - z.begin())*h;

    // График 1
    QCustomPlot* plot = ui->mPlot;
    plot->clearGraphs();

    plot->xAxis->setRange(-0.2, 1);
    plot->yAxis->setRange(-0.2, 1);

    plot->xAxis->setLabel("x");
    plot->yAxis->setLabel("u");
    plot->legend->setVisible(true);

    plot->addGraph();
    plot->graph(0)->addData(x, v);
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setName("Аналитическое решение v(x)");

    plot->addGraph();
    plot->graph(1)->addData(x, v2);
    plot->graph(1)->setPen(QPen(Qt::green));
    plot->graph(1)->setName("Аналитическое решение с половинным шагом v2(x)");

    plot->addGraph();
    plot->graph(2)->addData(x, z);
    plot->graph(2)->setPen(QPen(Qt::red));
    plot->graph(2)->setName("Модуль разности численных решений |v(x)-v2(x)|");

    plot->replot();

    if (ui->tButtonSolnInfo->isVisible())
        ui->tButtonSolnInfo->close();
    windowLabel->setText("Для решения задачи использована равномерная сетка с числом разбиейний n ="+ QString::number(n) +".\n"
                         "Задача должна быть решена с погрешностью не более Ɛ = 0.5*10^(-6).\n"
                         "Задача решена с погрешностью Ɛ₁ = " + QString::number(max_err) + ".\n"
                         "Максимальное разность численных решений в общих узлах сетки наблюдается в точке x = "+QString::number(max_err_pos)+".");
    ui->tButtonSolnInfo->click();
}

void MainWindow::on_tButtonSolnInfo_clicked()
{
    if (window != nullptr && !window->isHidden())
        window->hide();

    window->show();
}


void MainWindow::on_mButtonSolnInfo_pressed()
{
    if (window != nullptr && !window->isHidden())
        window->hide();

    window->show();
}

