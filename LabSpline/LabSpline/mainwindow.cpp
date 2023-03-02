#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cmath>
#include <limits>
#include <QMessageBox>
#include <QTableView>
#include "spline.h"

using std::sin;
using std::cos;
using std::abs;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    info = new QWidget();
    info->setWindowTitle("О программе");
    info->setGeometry(300, 300, 560, 370);
    solutionInfo = new QWidget();
    solutionInfo->setGeometry(300, 300, 350, 250);
    infoLabel = new QLabel(info);
    solutionInfoLabel = new QLabel(solutionInfo);
    solutionInfoLabel->setText("Построений ещё не было.\nВведите данные и нажмите \"Решить!\",\nЧтобы получить справку о решении.");
    infoLabel->setText("Программа выполнена в рамках лабораторной работы на тему\n\
\"Интерполяция кубическими сплайнами\". По заданным сетке и функции программа\n\
строит интерполяционный кубический сплайн на сетке.\n\
\n\
Есть 3 варианта задач в зависимости от интерполируемой функции F(x):\n\
1) F(x) = -|x|³ + 3x²;\n\
2) F(x) = sin(x) / (1 + x²);\n\
3) F(x) = sin(x) / (1 + x²) + cos(tx), t∈R - параметр.\n\
\n\
Пользователь вводит параметры:\n\
a, b - границы отрезка [a; b], на котором интерполируется F(x);\n\
S''(a), S''(b) - значения второй производной сплайна S(x) на концах отрезка [a; b];\n\
n - количество участков равномерной сетки, на которые делится отрезок [a; b];\n\
t - параметр в функции F(x) (доступен, если выбран 3 вариант задачи).\n\
\n\
После нажатия кнопки \"Решить!\" строятся графики F(x), S(x), а также их первой и второй\n\
производных, и в таблицы выводится значения этих функций на сетке и коэффициенты сплайнов.\n\
\n\
В таблице №1 приведены коэффициенты 'a_i', 'b_i', 'c_i' и 'd_i' для участка [x_(i-1); x_i],\n\
на котором S(x) имеет вид S_i(x) = a_i + b_i * (x - x_i) +c_i/2 * (x - x_i)² + d_i/6 * (x - x_i)³.\n\
В таблице №2 приведены значения F(x_j), S(x_j), F(x_j) - S(x_j), F'(x_j), S'(x_j), F'(x_j) - S'(x_j),\n\
F''(x_j), S''(x_j), F''(x_j) - S''(x_j).");
    func[0] = [](double x) {
        return -abs(x*x*x) + 3.0*x*x;
    };
    func[1] = [](double x) {
        return sin(x)/(1.0 + x*x);
    };
    func[2] = [=](double x) {
        return sin(x)/(1.0 + x*x) + cos(t*x);
    };
    funcD[0] = [](double x) {
        return -3.0*x*abs(x) + 6.0*x;
    };
    funcD[1] = [](double x) {
        double m = 1.0 + x*x;
        return (m*cos(x) - 2.0*x*sin(x))/(m*m);
    };
    funcD[2] = [=](double x) {
        double m = 1.0 + x*x;
        return (m*cos(x) - 2.0*x*sin(x))/(m*m) - t*sin(t*x);
    };
    funcD2[0] = [](double x) {
        return 6.0*(1.0 - abs(x));
    };
    funcD2[1] = [](double x) {
        double m = 1.0 + x*x;
        return -sin(x)/m - 4.0*x*cos(x)/(m*m) + (6.0*x*x - 2.0)*sin(x)/(m*m*m);
    };
    funcD2[2] = [=](double x) {
        double m = 1.0 + x*x;
        return -sin(x)/m - 4.0*x*cos(x)/(m*m) + (6.0*x*x - 2.0)*sin(x)/(m*m*m) - t*t*cos(t*x);
    };


    ui->table1->setModel(new QStandardItemModel(0,7,this));
    ui->table1->verticalHeader()->hide();
    ui->table1->setEditTriggers(QAbstractItemView::NoEditTriggers);
    QStandardItemModel *model = (QStandardItemModel*)ui->table1->model();
    model->setHeaderData(0, Qt::Horizontal, "i");
    model->setHeaderData(1, Qt::Horizontal, "x_(i-1)");
    model->setHeaderData(2, Qt::Horizontal, "x_i");
    model->setHeaderData(3, Qt::Horizontal, "a_i");
    model->setHeaderData(4, Qt::Horizontal, "b_i");
    model->setHeaderData(5, Qt::Horizontal, "c_i");
    model->setHeaderData(6, Qt::Horizontal, "d_i");

    ui->table2->setModel(new QStandardItemModel(0,11,this));
    ui->table2->verticalHeader()->hide();
    ui->table2->setEditTriggers(QAbstractItemView::NoEditTriggers);
    model = (QStandardItemModel*)ui->table2->model();
    model->setHeaderData(0, Qt::Horizontal, "j");
    model->setHeaderData(1, Qt::Horizontal, "x_j");
    model->setHeaderData(2, Qt::Horizontal, "F(x_j)");
    model->setHeaderData(3, Qt::Horizontal, "S(x_j)");
    model->setHeaderData(4, Qt::Horizontal, "F(x_j)-S(x_j)");
    model->setHeaderData(5, Qt::Horizontal, "F'(x_j)");
    model->setHeaderData(6, Qt::Horizontal, "S'(x_j)");
    model->setHeaderData(7, Qt::Horizontal, "F'(x_j)-S'(x_j)");
    model->setHeaderData(8, Qt::Horizontal, "F''(x_j)");
    model->setHeaderData(9, Qt::Horizontal, "S''(x_j)");
    model->setHeaderData(10, Qt::Horizontal, "F''(x_j)-S''(x_j)");

    ui->plot->setInteraction(QCP::iRangeDrag, true);
    ui->plot->setInteraction(QCP::iRangeZoom, true);
    ui->labelT->hide();
    ui->lineEditT->hide();
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_buttonSolve_clicked()
{
    int taskIdx = ui->comboBoxFunction->currentIndex();
    auto f = func[taskIdx];
    auto fD1 = funcD[taskIdx];
    auto fD2 = funcD2[taskIdx];
    bool ok;
    int n;
    double a,b,s2a,s2b;

    n = ui->lineEditN->text().toInt(&ok);
    if (ok == false || n <= 0) {
        QMessageBox::critical(this, "Ошибка!", "'n' должно быть натуральным числом");
        return;
    }
    a = ui->lineEditA->text().toDouble(&ok);
    if (ok == false) {
        QMessageBox::critical(this, "Ошибка!", "'a' должно быть действительным числом");
        return;
    }
    b = ui->lineEditB->text().toDouble(&ok);
    if (ok == false) {
        QMessageBox::critical(this, "Ошибка!", "'b' должно быть действительным числом");
        return;
    }
    if (a >= b) {
        QMessageBox::critical(this, "Ошибка!", "'a' должно быть меньше 'b'");
        return;
    }
    s2a = ui->lineEditSa->text().toDouble(&ok);
    if (ok == false) {
        QMessageBox::critical(this, "Ошибка!", "S''(a) должно быть действительным числом");
        return;
    }
    s2b = ui->lineEditSb->text().toDouble(&ok);
    if (ok == false) {
        QMessageBox::critical(this, "Ошибка!", "S''(b) должно быть действительным числом");
        return;
    }

    if (taskIdx == 2) {
        t = ui->lineEditT->text().toDouble(&ok);
        if (ok == false || t <= 0.0) {
            QMessageBox::critical(this, "Ошибка!", "'t' должно быть положительным числом");
            return;
        }
    }

    x.resize(n+1);
    double h = (b-a)/n;
    double xc = a;
    for (int i = 0; i <= n; ++i, xc += h)
        x[i] = xc;
    auto splines = findSplineCoefficients(x, f, s2a, s2b);

    // ЗАПОЛНЕНИЕ ДАННЫХ ДЛЯ ГРАФИКОВ
    cubicSpline spl = splines[0];
    QVector<double> plotX(50*n + 1),
                    plotS(50*n + 1), plotF(50*n + 1),
                    plotS1(50*n + 1), plotF1(50*n + 1),
                    plotS2(50*n + 1), plotF2(50*n + 1),
                    plotD(50*n + 1), plotD1(50*n + 1), plotD2(50*n + 1);
    int k = 1;
    double xcsq, xc1, xc2;
    double ymin = std::numeric_limits<double>::max(),
           ymax = std::numeric_limits<double>::min();
    for (int i = 0; i <= n; ++i) {
        ymin = std::min(ymin, f(x[i]));
        ymax = std::max(ymax, f(x[i]));
    }

    plotX[0] = a;
    plotS[0] = plotF[0] = f(x[0]);
    plotF1[0] = fD1(x[0]);
    plotS1[0] = spl.b - 2.0*spl.c*h + 3.0*spl.d*h*h;
    plotF2[0] = fD2(x[0]);
    plotS2[0] = 2.0*spl.c - 6.0*spl.d*h;
    plotD[0] = plotD1[0] = plotD2[0] = 0.0;
    h /= 50.0;
    xc1 = a + h;
    for (int i = 0; i < n; ++i) {
        spl = splines[i];
        xc2 = x[i] - x[i+1] + h;
        for (int j = 0; j < 50; ++j, ++k, xc1 += h, xc2 += h) {
            plotX[k] = xc1;
            xcsq = xc2 * xc2;
            plotF[k] = f(xc1);
            plotS[k] = spl.a + spl.b*xc2 + spl.c*xcsq + spl.d*xcsq*xc2;
            plotF1[k] = fD1(xc1);
            plotS1[k] = spl.b + 2.0*spl.c*xc2 + 3.0*spl.d*xcsq;
            plotF2[k] = fD2(xc1);
            plotS2[k] = 2.0*spl.c + 6.0*spl.d*xc2;
            plotD[k] = plotF[k] - plotS[k];
            plotD1[k] = plotF1[k] - plotS1[k];
            plotD2[k] = plotF2[k] - plotS2[k];
        }
    }

    // ПОСТРОЕНИЕ ГРАФИКОВ
    auto plot = ui->plot;
    plot->clearGraphs();
    for (int i = 0; i < 9; ++i)
        ui->plot->addGraph();
    plot->xAxis->setRange(a, b);
    plot->yAxis->setRange(ymin, ymax);
    plot->graph(0)->addData(plotX, plotF);
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setName("F(x)");
    plot->graph(1)->addData(plotX, plotS);
    plot->graph(1)->setPen(QPen(Qt::yellow));
    plot->graph(1)->setName("S(x)");
    plot->graph(2)->addData(plotX, plotF1);
    plot->graph(2)->setPen(QPen(Qt::red));
    plot->graph(2)->setName("F'(x)");
    plot->graph(3)->addData(plotX, plotS1);
    plot->graph(3)->setPen(QPen(Qt::green));
    plot->graph(3)->setName("S'(x)");
    plot->graph(4)->addData(plotX, plotF2);
    plot->graph(4)->setPen(QPen(Qt::magenta));
    plot->graph(4)->setName("F''(x)");
    plot->graph(5)->addData(plotX, plotS2);
    plot->graph(5)->setPen(QPen(Qt::gray));
    plot->graph(5)->setName("S''(x)");
    plot->graph(6)->addData(plotX, plotD);
    plot->graph(6)->setPen(QPen(Qt::darkCyan));
    plot->graph(6)->setName("F(x) - S(x)");
    plot->graph(7)->addData(plotX, plotD1);
    plot->graph(7)->setPen(QPen(Qt::black));
    plot->graph(7)->setName("F'(x) - S'(x)");
    plot->graph(8)->addData(plotX, plotD2);
    plot->graph(8)->setPen(QPen(Qt::darkGreen));
    plot->graph(8)->setName("F''(x) - S''(x)");

    if (!ui->checkBoxF->isChecked())
        plot->graph(0)->setVisible(false);
    if (!ui->checkBoxS->isChecked())
        plot->graph(1)->setVisible(false);
    if (!ui->checkBoxF1->isChecked())
        plot->graph(2)->setVisible(false);
    if (!ui->checkBoxS1->isChecked())
        plot->graph(3)->setVisible(false);
    if (!ui->checkBoxF2->isChecked())
        plot->graph(4)->setVisible(false);
    if (!ui->checkBoxS2->isChecked())
        plot->graph(5)->setVisible(false);
    if (!ui->checkBoxD->isChecked())
        plot->graph(6)->setVisible(false);
    if (!ui->checkBoxD1->isChecked())
        plot->graph(7)->setVisible(false);
    if (!ui->checkBoxD2->isChecked())
        plot->graph(8)->setVisible(false);

    plot->legend->setVisible(true);
    plot->replot();

    // ЗАПОЛНЕНИЕ ТАБЛИЦ
    auto model1 = (QStandardItemModel*)ui->table1->model();
    auto model2 = (QStandardItemModel*)ui->table2->model();
    model1->setRowCount(n);
    model2->setRowCount(50*n + 1);
    for (int i = 0; i < n; ++i) {
        spl = splines[i];
        model1->setData(model1->index(i,0), i+1);
        model1->setData(model1->index(i,1), x[i]);
        model1->setData(model1->index(i,2), x[i+1]);
        model1->setData(model1->index(i,3), spl.a);
        model1->setData(model1->index(i,4), spl.b);
        model1->setData(model1->index(i,5), spl.c);
        model1->setData(model1->index(i,6), spl.d);
    }
    for (int i = 0; i <= 50*n; ++i) {
        model2->setData(model2->index(i,0), i);
        model2->setData(model2->index(i,1), QString::number(plotX[i], 'f', 6));
        model2->setData(model2->index(i,2), QString::number(plotF[i], 'f', 8));
        model2->setData(model2->index(i,3), QString::number(plotS[i], 'f', 8));
        model2->setData(model2->index(i,4), QString::number(plotF[i] - plotS[i], 'e'));
        model2->setData(model2->index(i,5), QString::number(plotF1[i], 'f', 8));
        model2->setData(model2->index(i,6), QString::number(plotS1[i], 'f', 8));
        model2->setData(model2->index(i,7), QString::number(plotF1[i] - plotS1[i], 'e'));
        model2->setData(model2->index(i,8), QString::number(plotF2[i], 'f', 8));
        model2->setData(model2->index(i,9), QString::number(plotS2[i], 'f', 8));
        model2->setData(model2->index(i,10), QString::number(plotF2[i] - plotS2[i], 'e'));
    }
    ui->table1->resizeColumnsToContents();
    ui->table2->resizeColumnsToContents();

    // ЗАПОЛНЕНИЕ И ВЫВОД СПРАВКИ О РЕШЕНИИ
    QString solnInfoString;
    double maxD, maxD1, maxD2;
    int maxDIdx, maxD1Idx, maxD2Idx;
    maxD = maxD1 = maxD2 = 0.0;
    maxDIdx = maxD1Idx = maxD2Idx = 0;
    for (int i = 0; i <= 50*n; ++i) {
        if (abs(plotD[i]) > maxD) {
            maxD = abs(plotD[i]);
            maxDIdx = i;
        }
        if (abs(plotD1[i]) > maxD1) {
            maxD1 = abs(plotD1[i]);
            maxD1Idx = i;
        }
        if (abs(plotD2[i]) > maxD2) {
            maxD2 = abs(plotD2[i]);
            maxD2Idx = i;
        }
    }
    solnInfoString = "Справка\n\n\
Сетка сплайна: n = " + QString::number(n) +
"\nКонтрольная сетка: N = " + QString::number(50*n) +
"\n\nПогрешность сплайна на контрольной сетке:\n" +
"max {j=0..N} |F(x_j) - S(x_j)| = " + QString::number(maxD) +
" при x = " + QString::number(plotX[maxDIdx], 'f', 8) +
"\nПогрешность производной на контрольной сетке:\n" +
"max {j=0..N} |F'(x_j) - S'(x_j)| = " + QString::number(maxD1) +
" при x = " + QString::number(plotX[maxD1Idx], 'f', 8) +
"\nПогрешность второй производной на контрольной сетке:\n" +
"max {j=0..N} |F''(x_j) - S''(x_j)| = " + QString::number(maxD2) +
" при x = " + QString::number(plotX[maxD2Idx], 'f', 8);

    if (solutionInfo->isVisible())
        solutionInfo->hide();
    solutionInfoLabel->setText(std::move(solnInfoString));
    solutionInfo->show();
}

void MainWindow::on_comboBoxFunction_activated(int index)
{
    if (index == 2) {
        ui->lineEditT->show();
        ui->labelT->show();
    }
    else {
        ui->lineEditT->hide();
        ui->labelT->hide();
    }
}


void MainWindow::on_pushButton_2_clicked()
{
    if (!info->isHidden())
        info->hide();
    info->show();
}


void MainWindow::on_checkBoxF_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(0);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_checkBoxS_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(1);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_checkBoxF1_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(2);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_checkBoxS1_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(3);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}




void MainWindow::on_checkBoxF2_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(4);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_checkBoxS2_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(5);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}




void MainWindow::on_checkBoxD_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(6);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_checkBoxD1_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(7);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}

void MainWindow::on_checkBoxD2_clicked()
{
    if (ui->plot->graphCount()) {
        auto graph = ui->plot->graph(8);
        graph->setVisible(!graph->visible());
        ui->plot->replot();
    }
}


void MainWindow::on_buttonSolnInfo_clicked()
{
    if (solutionInfo->isVisible())
        solutionInfo->hide();
    solutionInfo->show();
}

