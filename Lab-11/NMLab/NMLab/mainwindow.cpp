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
    initializeBackground(plotWidget, backgroundWidget);
    plotWidget->setInteraction(QCP::iRangeZoom,true);
    plotWidget->setInteraction(QCP::iRangeDrag,true);
}

void MainWindow::initializeBackground(QWidget *widget, QWidget *backgroundWidget) {
    QPalette pal(palette());
    pal.setColor(QPalette::Window, Qt::gray);
    backgroundWidget->setAutoFillBackground(true);
    backgroundWidget->setPalette(pal);
    backgroundWidget->setGeometry(widget->x() - 1, widget->y() - 1, widget->width() + 2, widget->height() + 2);
}

void MainWindow::initializeMainTaskTable(QTableView *table, size_t rows)
{
    QStandardItemModel *model = new QStandardItemModel(rows,12,this);

    model->setHeaderData(0, Qt::Horizontal, "i");
    model->setHeaderData(1, Qt::Horizontal, "x_i");
    model->setHeaderData(2, Qt::Horizontal, "v_i");
    model->setHeaderData(3, Qt::Horizontal, "v_2_i");
    model->setHeaderData(4, Qt::Horizontal, "v_2_i - v_i");
    model->setHeaderData(5, Qt::Horizontal, "(v')_i");
    model->setHeaderData(6, Qt::Horizontal, "(v')_2_i");
    model->setHeaderData(7, Qt::Horizontal, "(v')_2_i - (v')_i");
    model->setHeaderData(8, Qt::Horizontal, "|S|");
    model->setHeaderData(9, Qt::Horizontal, "h_(i-1)");
    model->setHeaderData(10, Qt::Horizontal, "C+");
    model->setHeaderData(11, Qt::Horizontal, "C-");

    table->setModel(model);
    //table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
}

void MainWindow::on_progInfoButton_clicked()
{
    if (progInfo != nullptr && !progInfo->isHidden())
        progInfo->hide();

    progInfo->show();
}

void MainWindow::on_tableInfoButton_pressed()
{
    if (tableInfo != nullptr && !tableInfo->isHidden())
        tableInfo->hide();

    tableInfo->show();
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->setWindowTitle("Приложение по задаче 11. Шипилов Артём, 382003_3");

    initializeMainTaskTable(ui->main2Table, 0);
    ui->main2Table->setEditTriggers(QAbstractItemView::NoEditTriggers);

    progInfo    = new QWidget;
    tableInfo   = new QWidget;
    main2SolutionInfo = new QWidget;

    progInfoLabel    = new QLabel(progInfo);
    tableInfoLabel   = new QLabel(tableInfo);
    main2SolutionInfoLabel = new QLabel(main2SolutionInfo);

    progInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    tableInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    main2SolutionInfoLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

    initializeOutputInfoWidget(progInfo, progInfoLabel);
    initializeOutputInfoWidget(tableInfo, tableInfoLabel);
    initializeOutputInfoWidget(main2SolutionInfo, main2SolutionInfoLabel);

    progInfoLabel->setText(
"Это приложение выполнено в рамках задания 11 по численным методам.\n\
Суть задания - решить задачу Коши для ОДУ 2 порядка, используя метод Рунге-Кутты 4 порядка,\n\
показать зависимость решения от параметров модели, начальных условий и параметров численного метода,\n\
указать зависимость свойств решения, наиболее интересных с прикладной точки зрения, от параметров.\n\n\
Уравнение описывает движение математического маятника без трения в поле силы тяжести.\n\
'x'     - время [с].\n\
'u(x)'  - угол отклонения маятника по часовой стрелке от направленной вниз вертикали.\n\
'L'     - длина нити [м] (L > 0).\n\
'g'     - ускорение свободного падения [м/c^2] (принимается равным 9.8 м/c^2).\n\n\
Параметры задачи Коши:\n\
'u0'    - начальный угол отклонения от нисходящей вертикали [рад].\n\
'(u')0' - начальная угловая скорость маятника [рад/c].\n\n\
Параметры численного метода:\n\
'x_max' - правая граница отрезка численного решения (по времени).\n\
'ε_x'   - точность выхода на правую границу отрезка численного решения (по времени)\n\
(счёт прекращается, если для текущего времени 'x' верно: 'x > x_max - ε_x').\n\
'h'     - начальное значение шага интегрирования (по времени).\n\
'N_max' - максимальное число шагов метода (счёт прекращается при достижении шага с этим номером).\n\
'ε_ν'   - параметр контроля (участвует в оценке контрольного слагаемого и принятии решения об изменении шага).\n\n\
Критерии остановки счёта: выход на правую границу с заданной точностью, достижение заданного числа шагов метода.\n\
Метод контроля погрешности: счёт с половинным шагом и оценка величины '|S| = ||v_2_n - v_n||/15'.\n\
'|S| > ε_ν'            : полученные на этом шаге значения не принимаются; шаг уменьшается вдвое, .\n\
'ε_ν/32 <= |S| <= ε_ν': принимается значение, полученное с половинным шагом; шаг не меняется.\n\
'|S| < ε_ν/32'        : принимается значение, полученное с половинным шагом; шаг удваивается.\n\
");
    progInfo->setWindowTitle("О задаче");

    tableInfoLabel->setText("\
Будем далее полагать за 'v' численное решение, за ' v' ' - численно найденную производную численного решения.\n\n\
Смысл столбцов таблицы:\n\
1) 'i' - номер шага метода.\n\
2) 'x_i'      - значение 'x' текущей точки численного решения.\n\
3) 'v_i'      - значение 'v' в точке 'x_i', которое было принято методом.\n\
4) 'v_2_i'    - значение 'v' в точке 'x_i', полученное из точки 'x_(i-1)' двойным счётом с половинным шагом.\n\
    На 0-м шаге метода полагается равным 'u0' по определению.\n\
5) 'v_2_i     - v_i' - разность 'v_2_i' и 'v_i'.\n\
6) '(v')_i'   - значение ' v' ' в точке 'x_i', которое было принято методом.\n\
7) '(v')_2_i' - значение ' v' ' в точке 'x_i', полученное из точки 'x_(i-1)' двойным счётом с половинным шагом.\n\
  На 0-м шаге метода полагается равным '(u')0' по определению.\n\
8) '(v')_2_i - (v')_i' - разность '(v')_2_i' и '(v')_i'.\n\
9) '|S|'  - '|S| = ||v_i - v_2_i||/15' - величина, оцениваемая при принятии решения о продолжении счёта.\n\
10) 'h_(i-1)' - i >= 1: шаг численного интегрирования, 'h_(i-1) = x_i - x_(i-1)';\n\
                   i = 0: начальный шаг интегрирования.\n\
11) 'C+'  - количество удвоений шага интегрирования, выполненных на предыдущем шаге метода.\n\
12) 'C-' - количество делений вдвое шага интегрирования, выполненных на предыдущем шаге метода.\n\
");
    tableInfo->setWindowTitle("О столбцах таблицы");

    initializePlotWidget(ui->main2Plot1, ui->main2Plot1Background);
    initializePlotWidget(ui->main2Plot2, ui->main2Plot2Background);
    initializePlotWidget(ui->main2Plot3, ui->main2Plot3Background);
    initializeBackground(ui->main2Table, ui->main2TableBackground);
    ui->main2Table->verticalHeader()->setVisible(false);
}


MainWindow::~MainWindow()
{
    ui->main2Plot1->clearGraphs();
    ui->main2Plot2->clearGraphs();
    for (QCPCurve *curve: phase)
        delete curve;
    phase.clear();

    ui->main2Plot1->clearItems();
    ui->main2Plot2->clearItems();
    ui->main2Plot3->clearItems();

    delete main2SolutionInfo;
    delete progInfo,
    delete tableInfo;

    delete ui;
}

bool MainWindow::checkInput(const QString &str_x0, const QString &str_u0, const QString &str_x_max, const QString &str_h,
                            const QString &str_x_eps, const QString &str_n_max, const QString &str_v_eps,
                            double &x0, double &u0, double &x_max, double &h, double &x_eps, int &n_max, double &v_eps)
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
    else if (str_n_max.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'n_max' должно быть непустым");
    else if (str_v_eps.isEmpty())
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

    result = false;
    if (x_max < 0.0)
        QMessageBox::critical(this, "Ошибка!", "Правая граница отрезка 'x_max' должна быть не меньше чем левая - 'x0' = 0");
    else if (h <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Шаг интегрирования 'h' должен быть положительным");
    else if (x_eps <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Точность выхода на правую границу 'eps_x' должна быть положительной");
    else if (n_max <= 0)
        QMessageBox::critical(this, "Ошибка!", "Максимальное число шагов метода 'nmax' должно быть положительным");
    else if (v_eps <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Параметр контроля 'eps_v' должен быть положительным");
    else if (x_max < x0)
        QMessageBox::critical(this, "Ошибка!", "Правая граница 'x_max' быть не меньше 'x0'");
    else
        result = true;


    return result;
}

template <typename value_t, size_t dim>
void fillSolutionInfoLabel(QLabel *label, const pair<vector<Entry<value_t, dim>>, int> &soln, double x_max, QString &solutionInfoString)
{
    auto &solution = soln.first;
    int numOfIterations = soln.second;

    size_t n = solution.size();
    double h_min = solution[0].h, h_max = h_min;
    double s_max = 0.0, s_min = numeric_limits<double>::max();
    size_t h_min_idx = 0, h_max_idx = 0;
    size_t s_max_idx = 0, s_min_idx = 0;
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
        if (abs(solution[i].s) > s_max) {
            s_max = abs(solution[i].s);
            s_max_idx = i;
        }
        if (i && abs(solution[i].s) < s_min) {
            s_min = abs(solution[i].s);
            s_min_idx = i;
        }
    }

    solutionInfoString += "\n\nРезультат расчёта:";
    solutionInfoString += "\nЧисло шагов метода: N = " + QString::number(numOfIterations);
    solutionInfoString += "\nТочность выхода на правую границу: " + QString::number(x_max - solution.back().x);
    solutionInfoString += "\nx_N = " + QString::number(solution.back().x) + ", v_N = " + QString::number(solution.back().v[0]);
    solutionInfoString += "\n\nМаксимум модуля ОЛП: " + QString::number(s_max) + " при x = "
                                                      + QString::number(solution[s_max_idx].x);
    solutionInfoString += "\nМинимум модуля ОЛП: "  + QString::number(s_min) + " при x = "
                                                      + QString::number(solution[s_min_idx].x);

    solutionInfoString += "\n\nОбщее число удвоений шага: " + QString::number(plus_count);
    solutionInfoString += "\nОбщее число делений шага: " + QString::number(minus_count);

    solutionInfoString += "\n\nМинимальный шаг: " + QString::number(h_min) + " при x = " + QString::number(solution[h_min_idx].x);
    solutionInfoString += "\nМаксимальный шаг: " + QString::number(h_max) + " при x = " + QString::number(solution[h_max_idx].x);

    label->setText(solutionInfoString);
}

bool MainWindow::checkMain2Input(const QString & str_u0_dash, const QString & str_L, double &u0_dash, double &L)
{
    bool result = false;

    if (str_u0_dash.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле ' u'0 ' должно быть непустым");
    else if (str_L.isEmpty())
        QMessageBox::critical(this, "Ошибка!", "Поле 'L' должно быть непустым");
    else result = true;

    if (result == false)
        return false;
    u0_dash = str_u0_dash.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле ' u'0 ' имеет некорректное значение");
        return false;
    }
    L = str_L.toDouble(&result);
    if (result == false)
    {
        QMessageBox::critical(this, "Ошибка!", "Поле 'L' имеет некорректное значение");
        return false;
    }

    result = false;
    if (L <= 0.0)
        QMessageBox::critical(this, "Ошибка!", "Значение 'L' должно быть положительным");
    else
        result = true;

    return result;

}

void MainWindow::setGraphLabels(double L, double x0, double u0, double ud0) {
    QCPItemText *tLabel1 = new QCPItemText(ui->main2Plot1);
    tLabel1->position->setCoords(x0, u0);
    tLabel1->setText(QString("(%1, %2, %3, %4)").arg(L).arg(x0).arg(u0).arg(ud0));
    tLabel1->setPen(Qt::NoPen);
    tLabel1->setFont(QFont(font().family(), 12));
    tLabel1->setColor(Qt::darkMagenta);

    QCPItemText *tLabel2 = new QCPItemText(ui->main2Plot2);
    tLabel2->position->setCoords(x0, ud0);
    tLabel2->setText(QString("(%1, %2, %3, %4)").arg(L).arg(x0).arg(u0).arg(ud0));
    tLabel2->setPen(Qt::NoPen);
    tLabel2->setFont(QFont(font().family(), 12));
    tLabel2->setColor(Qt::darkMagenta);

    QCPItemText *tLabel3 = new QCPItemText(ui->main2Plot3);
    tLabel3->position->setCoords(u0, ud0);
    tLabel3->setText(QString("(%1, %2, %3, %4)").arg(L).arg(x0).arg(u0).arg(ud0));
    tLabel3->setPen(Qt::NoPen);
    tLabel3->setFont(QFont(font().family(), 12));
    tLabel3->setColor(Qt::darkMagenta);
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
            str_L     = ui->main2LineEdit_L->text(),
            str_n_max  = ui->main2LineEdit_nmax  ->text(),
            str_v_eps  = ui->main2LineEdit_v_eps ->text();

    double x0, u0, u0_dash, x_max, h, x_eps, v_eps, L;
    int n_max;
    if (!checkInput(str_x0, str_u0, str_x_max, str_h, str_x_eps, str_n_max, str_v_eps,
                        x0,     u0,     x_max,     h,     x_eps,     n_max,     v_eps)
        ||
        !checkMain2Input(str_u0_dash, str_L, u0_dash, L))
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
        return -(9.8*(1/L)*sin(y[0]));
    };
    auto soln = ivp_step_adjust(rk4_system2, Diff_equation<Point<2>,2>{f1,f2}, x0, Point<2>{u0,u0_dash/L}, h, x_max, x_eps, v_eps, n_max);
    auto &solution = soln.first;

// Заполнение таблицы
    size_t n = solution.size();
    auto model = (QStandardItemModel*)ui->main2Table->model();
    model->setRowCount(n);
    for (size_t i = 0; i < n; ++i) {
        model->setData(model->index(i,0), i);
        model->setData(model->index(i,1), QString::number(solution[i].x, 'f', 10));
        model->setData(model->index(i,2), QString::number(solution[i].v[0], 'f', 10));
        model->setData(model->index(i,3), QString::number(solution[i].v2[0], 'f', 10));
        model->setData(model->index(i,4), QString::number(solution[i].v2[0] - solution[i].v[0], 'f', 10));
        model->setData(model->index(i,5), QString::number(solution[i].v[1], 'f', 10));
        model->setData(model->index(i,6), QString::number(solution[i].v2[1], 'f', 10));
        model->setData(model->index(i,7), QString::number(solution[i].v2[1] - solution[i].v[1], 'f', 10));
        model->setData(model->index(i,8), solution[i].s);
        model->setData(model->index(i,9), solution[i].h);
        model->setData(model->index(i,10), solution[i].c_plus);
        model->setData(model->index(i,11), solution[i].c_minus);
    }
    ui->main2Table->resizeColumnsToContents();

// Построение графиков

    QVector<double> x(n), v(n), v2(n), v_dash(n), v2_dash(n);
    for (size_t i = 0; i < n; ++i)
    {
        x[i]        = solution[i].x;
        v[i]        = solution[i].v[0];
        v2[i]       = solution[i].v2[0];
        v_dash[i]   = solution[i].v[1];
        v2_dash[i]  = solution[i].v2[1];
    }

// График численного решения
    QCustomPlot *plot1 = ui->main2Plot1;

    double v_min = *min_element(v.begin(), v.end());

    plot1->xAxis->setRange(x0, x_max);
    plot1->yAxis->setRange(v_min, v_min + (x_max - x0));

    plot1->xAxis->setLabel("x");
    plot1->yAxis->setLabel("u");

    plot1->addGraph();
    plot1->graph(2*graphCount)->addData(x, v);
    plot1->graph(2*graphCount)->setPen(QPen(Qt::blue));

    plot1->addGraph();
    plot1->graph(2*graphCount + 1)->addData(x, v2);
    plot1->graph(2*graphCount + 1)->setPen(QPen(Qt::green));

// График производной численного решения
    QCustomPlot *plot2 = ui->main2Plot2;

    v_min = *min_element(v_dash.begin(), v_dash.end());

    plot2->xAxis->setRange(x0, x_max);
    plot2->yAxis->setRange(v_min, v_min + (x_max - x0));

    plot2->xAxis->setLabel("x");
    plot2->yAxis->setLabel("u'");

    plot2->addGraph();
    plot2->graph(2*graphCount)->addData(x, v_dash);
    plot2->graph(2*graphCount)->setPen(QPen(Qt::blue));

    plot2->addGraph();
    plot2->graph(2*graphCount + 1)->addData(x, v2_dash);
    plot2->graph(2*graphCount + 1)->setPen(QPen(Qt::green));

// Фазовый портрет
    QCustomPlot *plot3 = ui->main2Plot3;

    double xmin = *min_element(v.begin(), v.end()),
           xmax = *max_element(v.begin(), v.end());

    plot3->xAxis->setRange(u0 - (xmax - xmin)/2.0, u0 + (xmax - xmin)/2.0);
    plot3->yAxis->setRange(u0_dash - (xmax - xmin)/2.0, u0_dash + (xmax - xmin)/2.0);

    plot3->xAxis->setLabel("u");
    plot3->yAxis->setLabel("u'");

    phase.push_back(new QCPCurve(plot3->xAxis, plot3->yAxis));
    phase.back()->setData(v, v_dash);

    setGraphLabels(L, x0, u0, u0_dash);
    plot1->replot();
    plot2->replot();
    plot3->replot();
    ++graphCount;

// Создание справки и её отображение
    main2IsSolvedOnce = true;

    QString solutionInfoString;
    solutionInfoString += "№ варианта задания: 1";
    solutionInfoString += "\nТип задачи: основная";
    solutionInfoString += "\nМетод Рунге-Кутта порядка 'p = 4'";
    solutionInfoString += "\nСпособы счёта - выход на правую границу ИЛИ выполнение заданного числа шагов метода";
    solutionInfoString += "\n\nx0 = " + QString::number(x0) + ", u0 = " + QString::number(u0);
    solutionInfoString += "\nx_max = " + QString::number(x_max) + ", ε_x = " + QString::number(x_eps);
    solutionInfoString += "\nh0 = " + QString::number(h) + ", N_max = " + QString::number(n_max);
    solutionInfoString += "\nКонтроль локальной погрешности включён";
    solutionInfoString += "\nε = " + QString::number(v_eps) + ", ε_min = " + QString::number(v_eps/32.0);
    fillSolutionInfoLabel(main2SolutionInfoLabel, soln, x_max, solutionInfoString);
    ui->main2ShowSolnInfoButton->click();
}


void MainWindow::on_main2ShowSolnInfoButton_clicked()
{
    if (!main2IsSolvedOnce)
        return;

    if (main2SolutionInfo != nullptr && !main2SolutionInfo->isHidden())
        main2SolutionInfo->hide();

    main2SolutionInfo->show();
}


void MainWindow::on_clearButton_pressed()
{
    ui->main2Plot1->clearGraphs();
    ui->main2Plot2->clearGraphs();
    for (QCPCurve *curve: phase)
        delete curve;
    phase.clear();

    ui->main2Plot1->clearItems();
    ui->main2Plot2->clearItems();
    ui->main2Plot3->clearItems();
    ui->main2Plot1->replot();
    ui->main2Plot2->replot();
    ui->main2Plot3->replot();
    graphCount = 0;
}

