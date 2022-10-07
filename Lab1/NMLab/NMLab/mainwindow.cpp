#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QDebug>
#include <Methods.h>
#include <qcustomplot.h>
#include <limits>
#include <QtMath>
using namespace std;

void MainWindow::initializeOutputInfoWidget(QWidget *parent, QLabel *label)
{
    parent->setWindowTitle("Выходные данные");
    parent->setGeometry(300, 300, 300, 300);
    QHBoxLayout *layout = new QHBoxLayout(parent);
    layout->setAlignment(parent, Qt::AlignCenter);
    layout->addWidget(label);
    parent->setLayout(layout);
}

void MainWindow::initializePlotWidget(QCustomPlot *plotWidget, QWidget *backgroundWidget)
{
    QPalette pal(palette());
    pal.setColor(QPalette::Background, Qt::gray);
    backgroundWidget->setAutoFillBackground(true);
    backgroundWidget->setPalette(pal);
    backgroundWidget->setGeometry(plotWidget->x() - 1, plotWidget->y() - 1, plotWidget->width() + 2, plotWidget->height() + 2);
    plotWidget->setInteraction(QCP::iRangeZoom,true);
    plotWidget->setInteraction(QCP::iRangeDrag,true);
}

void MainWindow::initializeTestTaskTable(QTableView *table, size_t rows)
{
    QStandardItemModel *model = new QStandardItemModel(rows,10,this);

    model->setHeaderData(0, Qt::Horizontal, "x_i");
    model->setHeaderData(1, Qt::Horizontal, "v_i");
    model->setHeaderData(2, Qt::Horizontal, "v_2i");
    model->setHeaderData(3, Qt::Horizontal, "v_i - v_2i");
    model->setHeaderData(4, Qt::Horizontal, "ОЛП");
    model->setHeaderData(5, Qt::Horizontal, "u_i");
    model->setHeaderData(6, Qt::Horizontal, "|u_i - v_i|");
    model->setHeaderData(7, Qt::Horizontal, "h_i");
    model->setHeaderData(8, Qt::Horizontal, "C+");
    model->setHeaderData(9, Qt::Horizontal, "C-");

    table->setModel(model);
}

void MainWindow::initializeMainTaskTable(QTableView *table, size_t rows)
{
    QStandardItemModel *model = new QStandardItemModel(rows,8,this);

    model->setHeaderData(0, Qt::Horizontal, "x_i");
    model->setHeaderData(1, Qt::Horizontal, "v_i");
    model->setHeaderData(2, Qt::Horizontal, "v_2i");
    model->setHeaderData(3, Qt::Horizontal, "v_i - v_2i");
    model->setHeaderData(4, Qt::Horizontal, "ОЛП");
    model->setHeaderData(5, Qt::Horizontal, "h_i");
    model->setHeaderData(6, Qt::Horizontal, "C+");
    model->setHeaderData(7, Qt::Horizontal, "C-");

    table->setModel(model);
}

void MainWindow::on_introTextButton_clicked()
{
    if (introReference != nullptr && !introReference->isHidden())
        introReference->hide();

    introReference->show();
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->setWindowTitle("Numerical Analysis Laboratory Work");

    initializeTestTaskTable(ui->testTable, 0);
    ui->testTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    initializeMainTaskTable(ui->main1Table, 0);
    ui->main1Table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    initializeMainTaskTable(ui->main2Table, 0);
    ui->main2Table->setEditTriggers(QAbstractItemView::NoEditTriggers);

    introReference    = new QWidget;
    testSolutionInfo  = new QWidget;
    main1SolutionInfo = new QWidget;
    main2SolutionInfo = new QWidget;

    introReferenceLabel    = new QLabel(introReference);
    testSolutionInfoLabel  = new QLabel(testSolutionInfo);
    main1SolutionInfoLabel = new QLabel(main1SolutionInfo);
    main2SolutionInfoLabel = new QLabel(main2SolutionInfo);

    testSolutionInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    main1SolutionInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    main2SolutionInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

    initializeOutputInfoWidget(introReference, introReferenceLabel);
    initializeOutputInfoWidget(testSolutionInfo , testSolutionInfoLabel);
    initializeOutputInfoWidget(main1SolutionInfo, main1SolutionInfoLabel);
    initializeOutputInfoWidget(main2SolutionInfo, main2SolutionInfoLabel);

    introReferenceLabel->setText(
    "Это приложение выполнено в рамках лабораторной работы №1 по численным методам.\n\
Оно содержит 3 вкладки, каждая из них посвящена численному решению задачи Коши для конкретного ОДУ.\n\n\
Пользователь может задать начальные условия (точку и значение функции в ней) и параметры численного решения\n\
(правая граница отрезка 'x_max', точность выхода на правую границу 'eps_x', начальная величина шага 'h', \n\
для методов с регулировкой шага - максимальное число итераций 'N_max' и параметр контроля 'ε').\n\n\
После введения всех значений и нажатия кнопки 'Решить!' заполняется таблица с данными о решении и строится \n\
график численного решения (в тестовой задаче - ещё и истинное решение).\n\n\
Если введённые значения некорректны, выведется окно с сообщением об ошибке.");
    introReference->setWindowTitle("О программе");

    initializePlotWidget(ui->testPlot,  ui->testPlotBackground);
    initializePlotWidget(ui->main1Plot, ui->main1PlotBackground);
    initializePlotWidget(ui->main2Plot1, ui->main2Plot1Background);
    initializePlotWidget(ui->main2Plot2, ui->main2Plot2Background);
    initializePlotWidget(ui->main2Plot3, ui->main2Plot3Background);
    phase = new QCPCurve(ui->main2Plot3->xAxis, ui->main2Plot3->yAxis);
    phase->setName("Фазовая траектория");

    ui->testRadioButton ->click();
    ui->main1RadioButton->click();
    ui->main2RadioButton->click();

    ui->tabWidget->setCurrentIndex(0);
    ui->introTextButton->click();
    introReference->activateWindow();
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_testRadioButton_clicked()
{
    ui->testLineEdit_nmax->hide();
    ui->testPicture_nmax ->hide();

    ui->testLineEdit_v_eps ->hide();
    ui->testPicture_v_eps  ->hide();
}

void MainWindow::on_testRadioButton2_clicked()
{
    ui->testLineEdit_nmax->show();
    ui->testPicture_nmax ->show();

    ui->testLineEdit_v_eps ->show();
    ui->testPicture_v_eps  ->show();
}

bool MainWindow::checkInput(const QString &str_x0, const QString &str_u0, const QString &str_x_max, const QString &str_h,
                            const QString &str_x_eps, const QString &str_n_max, const QString &str_v_eps,
                            double &x0, double &u0, double &x_max, double &h, double &x_eps, int &n_max, double &v_eps,
                            QRadioButton *&radioButton)
{

    bool result = false;

    if (str_x0.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'x0' должно быть непустым");
    else if (str_u0.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'u0' должно быть непустым");
    else if (str_x_max.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'x_max' должно быть непустым");
    else if (str_h.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'h' должно быть непустым");
    else if (str_x_eps.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'eps_x' должно быть непустым");
    else if (radioButton->isChecked() && str_n_max.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'n_max' должно быть непустым");
    else if (radioButton->isChecked() && str_v_eps.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'eps_v' должно быть непустым");
    else
        result = true;

    if (result == false)
        return false;
    x0    = str_x0.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'x0' имеет некорректное значение");
        return false;
    }
    u0    = str_u0.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'u0' имеет некорректное значение");
        return false;
    }
    x_max = str_x_max.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'x_max' имеет некорректное значение");
        return false;
    }
    h     = str_h.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'h' имеет некорректное значение");
        return false;
    }
    x_eps = str_x_eps.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'eps_x' имеет некорректное значение");
        return false;
    }
    if(radioButton->isChecked())
    {
        n_max = str_n_max.toInt(&result);
        if (result == false)
        {
            QMessageBox::critical(this, "Ошибка!", "Поле 'n_max' имеет некорректное значение");
            return false;
        }
        v_eps = str_v_eps.toDouble(&result);
        if (result == false)
        {
            QMessageBox::critical(this, "Ошибка!", "Поле 'eps_v' имеет некорректное значение");
            return false;
        }
    }

    result = false;
    if (x_max < 0.0)
        QMessageBox::critical(this, "Ошибка!", "Правая граница отрезка 'x_max' должна быть не меньше чем левая - 'x0' = 0");
    else if (h <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Шаг интегрирования 'h' должен быть положительным");
    else if (x_eps <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Точность выхода на правую границу 'eps_x' должна быть положительной");
    else if (radioButton->isChecked() && n_max <= 0)
        QMessageBox::critical(this, "Ошибка!", "Максимальное число шагов метода 'nmax' должно быть положительным");
    else if (radioButton->isChecked() && v_eps <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Параметр контроля 'eps_v' должен быть положительным");
    else
        result = true;


    return result;
}

template <typename value_t, size_t dim>
void fillSolutionInfoLabel(QLabel *label, const vector<Entry<value_t, dim>> &solution, double x_max)
{
    size_t n = solution.size();
    double h_min = solution[0].h, h_max = h_min;
    double s_star_max = 0.0;
    size_t h_min_idx = 0, h_max_idx = 0;
    size_t plus_count = 0, minus_count = 0;

    for (size_t i = 0; i < n; ++i)
    {
        plus_count += solution[i].c_plus;
        minus_count += solution[i].c_minus;
        if (solution[i].h < h_min)
        {
            h_min_idx = i;
            h_min = solution[i].h;
        }
        if (solution[i].h > h_max)
        {
            h_max_idx = i;
            h_max = solution[i].h;
        }
        if (abs(solution[i].s_star) > s_star_max)
            s_star_max = abs(solution[i].s_star);
    }

    QString solutionInfoString;
    solutionInfoString += "Число шагов метода (N): " + QString::number(n-1);
    solutionInfoString += "\nx_max  -  x_N: " + QString::number(x_max - solution.back().x);
    solutionInfoString += "\n\nМаксимум модуля ОЛП: " + QString::number(s_star_max);

    solutionInfoString += "\n\nОбщее число удвоений шага: " + QString::number(plus_count);
    solutionInfoString += "\nОбщее число делений шага: " + QString::number(minus_count);

    solutionInfoString += "\n\nМинимальный шаг: " + QString::number(h_min) + " при x = " + QString::number(solution[h_min_idx].x);
    solutionInfoString += "\nМаксимальный шаг: " + QString::number(h_max) + " при x = " + QString::number(solution[h_max_idx].x);

    label->setText(solutionInfoString);
}


void MainWindow::on_testSolveButton_pressed()
{

// Ввод данных и проверка их корректности
    QString str_x0    = ui->testLineEdit_x0   ->text(),
            str_u0    = ui->testLineEdit_u0   ->text(),
            str_x_max = ui->testLineEdit_xmax ->text(),
            str_h     = ui->testLineEdit_h    ->text(),
            str_x_eps = ui->testLineEdit_x_eps->text(),
            str_n_max,
            str_v_eps;
    if(ui->testRadioButton2->isChecked())
    {
        str_n_max  = ui->testLineEdit_nmax  ->text();
        str_v_eps  = ui->testLineEdit_v_eps ->text();
    }

    double x0, u0, x_max, h, x_eps, v_eps;
    int n_max;
    if (!checkInput(str_x0, str_u0, str_x_max, str_h, str_x_eps, str_n_max, str_v_eps,
                        x0,     u0,     x_max,     h,     x_eps,     n_max,     v_eps, ui->testRadioButton2))
    {
        return;
    }


// Численное решение
    auto f = [](double x, double u)
    {
        return -2.5*u;
    };
    auto solution = (ui->testRadioButton2->isChecked())
            ?    ivp_step_adjust(rk4, Diff_equation<double,1>{f}, x0, u0, h, x_max, x_eps, v_eps, n_max)
            : ivp_no_step_adjust(rk4, Diff_equation<double,1>{f}, x0, u0, h, x_max, x_eps);

    size_t n = solution.size();
    QVector<double> u(n);
    for (size_t i = 0; i < n; ++i)
        u[i] = u0 * exp(-2.5 * (solution[i].x - x0));

    QVector<double> x(n), v(n), v2(n);
    for (size_t i = 0; i < n; ++i)
    {
        x[i] = solution[i].x;
        v[i] = solution[i].v;
        v2[i] = solution[i].v2;
    }


// Заполнение таблицы
    auto model = (QStandardItemModel*)ui->testTable->model();
    model->setRowCount(n);
    for (size_t i = 0; i < n; ++i)
    {
        model->setData(model->index(i,0), solution[i].x);
        model->setData(model->index(i,1), solution[i].v);
        model->setData(model->index(i,2), solution[i].v2);
        model->setData(model->index(i,3), solution[i].v - solution[i].v2);
        model->setData(model->index(i,4), solution[i].s_star);
        model->setData(model->index(i,5), u[i]);
        model->setData(model->index(i,6), abs(solution[i].v - u[i]));
        model->setData(model->index(i,7), solution[i].h);
        model->setData(model->index(i,8), solution[i].c_plus);
        model->setData(model->index(i,9), solution[i].c_minus);
    }


// Построение графиков
    QCustomPlot *plot = ui->testPlot;

    plot->clearGraphs();

    plot->xAxis->setRange(x0,x_max);
    plot->yAxis->setRange(*min_element(v.begin(),v.end()), *max_element(v.begin(), v.end()));

    plot->legend->setVisible(true);
    plot->xAxis->setLabel("x");
    plot->yAxis->setLabel("u");

    plot->addGraph();
    plot->graph(0)->addData(x, v);
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setName("Численное решение v(x)");

    plot->addGraph();
    plot->graph(1)->addData(x, v2);
    plot->graph(1)->setPen(QPen(Qt::green));
    plot->graph(1)->setName("Численное решение c половинным шагом v2(x)");

    plot->addGraph();
    plot->graph(2)->addData(x, u);
    plot->graph(2)->setPen(QPen(Qt::red));
    plot->graph(2)->setName("Истинное решение  u(x)");

    plot->replot();


// Создание справки и её отображение
    testIsSolvedOnce = true;

    QLabel *label = testSolutionInfoLabel;
    fillSolutionInfoLabel(label, solution, x_max);

    double err_max = abs(v[0] - u[0]);
    size_t err_max_idx = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (abs(v[i] - u[i]) > err_max)
        {
            err_max = abs(v[i] - u[i]);
            err_max_idx = i;
        }
    }
    label->setText(label->text() + "\n\nМаксимальная глобальная погрешность: " + QString::number(err_max) +
                                   " при x = " + QString::number(x[err_max_idx]));

    ui->testShowSolnInfoButton->click();
}

void MainWindow::on_testShowSolnInfoButton_clicked()
{
    if (!testIsSolvedOnce)
        return;

    if (testSolutionInfo != nullptr && !testSolutionInfo->isHidden())
        testSolutionInfo->hide();

    testSolutionInfo->show();
}


// ОСНОВНОЕ ЗАДАНИЕ 1

void MainWindow::on_main1SolveButton_clicked()
{
// Ввод данных и проверка их корректности
    QString str_x0    = ui->main1LineEdit_x0   ->text(),
            str_u0    = ui->main1LineEdit_u0   ->text(),
            str_x_max = ui->main1LineEdit_xmax ->text(),
            str_h     = ui->main1LineEdit_h    ->text(),
            str_x_eps = ui->main1LineEdit_x_eps->text(),
            str_n_max,
            str_v_eps;
    if(ui->main1RadioButton2->isChecked())
    {
        str_n_max  = ui->main1LineEdit_nmax  ->text();
        str_v_eps  = ui->main1LineEdit_v_eps ->text();
    }

    double x0, u0, x_max, h, x_eps, v_eps;
    int n_max;
    if (!checkInput(str_x0, str_u0, str_x_max, str_h, str_x_eps, str_n_max, str_v_eps,
                        x0,     u0,     x_max,     h,     x_eps,     n_max,     v_eps, ui->main1RadioButton2))
    {
        return;
    }


// Численное решение
    auto f = [](double x, double u)
    {
        return u*(log(x+1)/(x*x + 1.0)*u - 1.0 + u*u*sin(10.0*x));
    };
    auto solution = (ui->main1RadioButton2->isChecked())
            ?    ivp_step_adjust(rk4, Diff_equation<double,1>{f}, x0, u0, h, x_max, x_eps, v_eps, n_max)
            : ivp_no_step_adjust(rk4, Diff_equation<double,1>{f}, x0, u0, h, x_max, x_eps);


// Заполнение таблицы
    size_t n = solution.size();
    QVector<double> x(n), v(n), v2(n);
    for (size_t i = 0; i < n; ++i)
    {
        x[i] = solution[i].x;
        v[i] = solution[i].v;
        v2[i] = solution[i].v2;
    }

    auto model = (QStandardItemModel*)ui->main1Table->model();
    model->setRowCount(n);
    for (size_t i = 0; i < n; ++i) {
        model->setData(model->index(i,0), solution[i].x);
        model->setData(model->index(i,1), solution[i].v);
        model->setData(model->index(i,2), solution[i].v2);
        model->setData(model->index(i,3), solution[i].v - solution[i].v2);
        model->setData(model->index(i,4), solution[i].s_star);
        model->setData(model->index(i,5), solution[i].h);
        model->setData(model->index(i,6), solution[i].c_plus);
        model->setData(model->index(i,7), solution[i].c_minus);
    }


// Построение графиков
    QCustomPlot *plot = ui->main1Plot;
    plot->clearGraphs();

    double y_min = *min_element(v.begin(), v.end()),
           y_max = *max_element(v.begin(), v.end());

    plot->xAxis->setRange(0,x_max);
    plot->yAxis->setRange(y_min, y_max);

    plot->legend->setVisible(true);
    plot->xAxis->setLabel("x");
    plot->yAxis->setLabel("u");

    plot->addGraph();
    plot->graph(0)->addData(x, v);
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setName("Численное решение v(x)");

    plot->addGraph();
    plot->graph(1)->addData(x, v2);
    plot->graph(1)->setPen(QPen(Qt::green));
    plot->graph(1)->setName("Численное решение с половинным шагом v2(x)");

    plot->replot();


// Создание справки и её отображение
    main1IsSolvedOnce = true;

    fillSolutionInfoLabel(main1SolutionInfoLabel, solution, x_max);
    ui->main1ShowSolnInfoButton->click();

}

void MainWindow::on_main1RadioButton_clicked()
{
    ui->main1Picture_nmax  ->hide();
    ui->main1LineEdit_nmax ->hide();

    ui->main1Picture_v_eps ->hide();
    ui->main1LineEdit_v_eps->hide();
}

void MainWindow::on_main1RadioButton2_clicked()
{
    ui->main1Picture_nmax  ->show();
    ui->main1LineEdit_nmax ->show();

    ui->main1Picture_v_eps ->show();
    ui->main1LineEdit_v_eps->show();
}

void MainWindow::on_main1ShowSolnInfoButton_clicked()
{
    if (!main1IsSolvedOnce)
        return;

    if (main1SolutionInfo != nullptr && !main1SolutionInfo->isHidden())
        main1SolutionInfo->hide();

    main1SolutionInfo->show();
}

bool MainWindow::checkMain2Input(const QString & str_u0_dash, const QString & str_a, const QString & str_b, double &u0_dash, double &a, double &b)
{
    bool result = false;

    if (str_u0_dash.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле ' u'0 ' должно быть непустым");
    else if (str_a.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'a' должно быть непустым");
    else if (str_b.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'b' должно быть непустым");
    else result = true;

    if (result == false)
        return false;
    u0_dash = str_u0_dash.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле ' u'0 ' имеет некорректное значение");
        return false;
    }
    a = str_a.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'a' имеет некорректное значение");
        return false;
    }
    b = str_b.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'b' имеет некорректное значение");
        return false;
    }

    result = false;
    if (a <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Значение 'a' должно быть положительным");
    else if (b <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Значение 'b' должно быть положительным");
    else
        result = true;

    return result;

}

void MainWindow::on_main2SolveButton_clicked()
{
// Ввод данных и проверка их корректности
    QString str_x0    = ui->main2LineEdit_x0   ->text(),
            str_u0    = ui->main2LineEdit_u0   ->text(),
          str_u0_dash = ui->main2LineEdit_u0_dash->text(),
            str_x_max = ui->main2LineEdit_xmax ->text(),
            str_h     = ui->main2LineEdit_h    ->text(),
            str_x_eps = ui->main2LineEdit_x_eps->text(),
            str_a     = ui->main2LineEdit_a    ->text(),
            str_b     = ui->main2LineEdit_b    ->text(),
            str_n_max,
            str_v_eps;
    if(ui->main2RadioButton2->isChecked())
    {
        str_n_max  = ui->main2LineEdit_nmax  ->text();
        str_v_eps  = ui->main2LineEdit_v_eps ->text();
    }

    double x0, u0, u0_dash, x_max, h, x_eps, v_eps, a, b;
    int n_max;
    if (!checkInput(str_x0, str_u0, str_x_max, str_h, str_x_eps, str_n_max, str_v_eps,
                        x0,     u0,     x_max,     h,     x_eps,     n_max,     v_eps, ui->main2RadioButton2)
        ||
        !checkMain2Input(str_u0_dash, str_a, str_b, u0_dash, a, b))
    {
        return;
    }

// Численное решение
    auto f1 = [](double x, const Point<2> &y)
    {
        return y[1];
    };
    auto f2 = [&](double x, const Point<2> &y)
    {
        return -(a*y[1]*y[1] + b*sin(y[0]));
    };
    auto solution = (ui->main2RadioButton2->isChecked())
            ?    ivp_step_adjust(rk4_system2, Diff_equation<Point<2>,2>{f1,f2}, x0, Point<2>{u0,u0_dash}, h, x_max, x_eps, v_eps, n_max)
            : ivp_no_step_adjust(rk4_system2, Diff_equation<Point<2>,2>{f1,f2}, x0, Point<2>{u0,u0_dash}, h, x_max, x_eps);

// Заполнение таблицы
    size_t n = solution.size();
    auto model = (QStandardItemModel*)ui->main2Table->model();
    model->setRowCount(n);
    for (size_t i = 0; i < n; ++i) {
        model->setData(model->index(i,0), solution[i].x);
        model->setData(model->index(i,1), solution[i].v[0]);
        model->setData(model->index(i,2), solution[i].v2[0]);
        model->setData(model->index(i,3), solution[i].v[0] - solution[i].v2[0]);
        model->setData(model->index(i,4), solution[i].s_star);
        model->setData(model->index(i,5), solution[i].h);
        model->setData(model->index(i,6), solution[i].c_plus);
        model->setData(model->index(i,7), solution[i].c_minus);
    }
    ui->main2Table->setModel(model);

// Построение графиков

    QVector<double> x(n), v(n), v2(n), v_dash(n), v2_dash(n);
    for (size_t i = 0; i < n; ++i)
    {
        x[i]        = solution[i].x;
        v[i]        = solution[i].v[0];
        v2[i]      = solution[i].v2[0];
        v_dash[i]   = solution[i].v[1];
        v2_dash[i] = solution[i].v2[1];
    }

// График численного решения
    QCustomPlot *plot1 = ui->main2Plot1;
    plot1->clearGraphs();

    double y_min = *min_element(v.begin(), v.end()),
           y_max = *max_element(v.begin(), v.end());

    plot1->xAxis->setRange(0,x_max);
    plot1->yAxis->setRange(y_min, y_max);

    plot1->legend->setVisible(true);
    plot1->xAxis->setLabel("x");
    plot1->yAxis->setLabel("u");

    plot1->addGraph();
    plot1->graph(0)->addData(x, v);
    plot1->graph(0)->setPen(QPen(Qt::blue));
    plot1->graph(0)->setName("Численное решение v(x)");

    plot1->addGraph();
    plot1->graph(1)->addData(x, v2);
    plot1->graph(1)->setPen(QPen(Qt::green));
    plot1->graph(1)->setName("Численное решение с половинным шагом v2(x)");

    plot1->replot();

// График производной численного решения
    QCustomPlot *plot2 = ui->main2Plot2;
    plot2->clearGraphs();

    y_min = *min_element(v_dash.begin(), v_dash.end()),
    y_max = *max_element(v_dash.begin(), v_dash.end());

    plot2->xAxis->setRange(0,x_max);
    plot2->yAxis->setRange(y_min, y_max);

    plot2->legend->setVisible(true);
    plot2->xAxis->setLabel("x");
    plot2->yAxis->setLabel("u'");

    plot2->addGraph();
    plot2->graph(0)->addData(x, v_dash);
    plot2->graph(0)->setPen(QPen(Qt::blue));
    plot2->graph(0)->setName("Производная численного решения v(x)");

    plot2->addGraph();
    plot2->graph(1)->addData(x, v2_dash);
    plot2->graph(1)->setPen(QPen(Qt::green));
    plot2->graph(1)->setName("Производная численного решения с половинным шагом v2(x)");

    plot2->replot();


// Фазовый портрет
    QCustomPlot *plot3 = ui->main2Plot3;
    plot3->clearGraphs();

    double xmin = *min_element(v.begin(), v.end()),
           xmax = *max_element(v.begin(), v.end());

    plot3->xAxis->setRange(u0 - (xmax - xmin)/2.0, u0 + (xmax - xmin)/2.0);
    plot3->yAxis->setRange(u0_dash - (xmax - xmin)/2.0, u0_dash + (xmax - xmin)/2.0);

    plot3->legend->setVisible(true);
    plot3->xAxis->setLabel("u");
    plot3->yAxis->setLabel("u'");

    phase->setData(v, v_dash);

    plot3->replot();


// Создание справки и её отображение
    main2IsSolvedOnce = true;

    fillSolutionInfoLabel(main2SolutionInfoLabel, solution, x_max);
    ui->main2ShowSolnInfoButton->click();
}

void MainWindow::on_main2RadioButton_clicked()
{
    ui->main2Picture_nmax  ->hide();
    ui->main2LineEdit_nmax ->hide();

    ui->main2Picture_v_eps ->hide();
    ui->main2LineEdit_v_eps->hide();
}

void MainWindow::on_main2RadioButton2_clicked()
{
    ui->main2Picture_nmax  ->show();
    ui->main2LineEdit_nmax ->show();

    ui->main2Picture_v_eps ->show();
    ui->main2LineEdit_v_eps->show();
}

void MainWindow::on_main2ShowSolnInfoButton_clicked()
{
    if (!main2IsSolvedOnce)
        return;

    if (main2SolutionInfo != nullptr && !main2SolutionInfo->isHidden())
        main2SolutionInfo->hide();

    main2SolutionInfo->show();
}
