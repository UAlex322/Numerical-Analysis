#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QMessageBox>
#include <QStandardItemModel>
#include <QThread>
#include <string>
#include <cmath>
using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    this->setWindowTitle("Устинов Александр, 382003_3");
    graph = new Q3DSurface;
    plot = QWidget::createWindowContainer(graph);
    ui->horizontalLayout->addWidget(plot, 1);
    ui->verticalLayout->setAlignment(Qt::AlignTop);
    graph->axisX()->setTitle("x");
    graph->axisY()->setTitle("u(x,y)");
    graph->axisZ()->setTitle("y");
    graph->axisX()->setTitleVisible(true);
    graph->axisY()->setTitleVisible(true);
    graph->axisZ()->setTitleVisible(true);

    solutionInfo = new QWidget;
    solutionInfo->setWindowTitle("Справка о последнем решении");
    solutionInfo->setGeometry(200, 200, 520, 400);
    solutionLabel = new QLabel(solutionInfo);
    solutionLabel->setGeometry(10, 10, 500, 395);
    // solutionLabel->setAlignment(Qt::AlignHCenter);
    solutionLabel->setText("Итерационный метод ещё ни разу не был запущен.");

    uSeries  = new QSurface3DSeries;
    v0Series = new QSurface3DSeries;
    vnSeries = new QSurface3DSeries;
    znSeries = new QSurface3DSeries;
    uSeries->setBaseColor(Qt::green);
    uSeries->setDrawMode(QSurface3DSeries::DrawSurface);
    uSeries->setFlatShadingEnabled(true);
    v0Series->setBaseColor(Qt::yellow);
    v0Series->setDrawMode(QSurface3DSeries::DrawSurface);
    v0Series->setFlatShadingEnabled(true);
    vnSeries->setBaseColor(Qt::blue);
    vnSeries->setDrawMode(QSurface3DSeries::DrawSurface);
    vnSeries->setFlatShadingEnabled(true);
    znSeries->setBaseColor(Qt::red);
    znSeries->setDrawMode(QSurface3DSeries::DrawSurface);
    znSeries->setFlatShadingEnabled(true);

    solutionLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);

    ui->tableViewU->setModel(new QStandardItemModel);
    ui->tableViewU->setFont(QFont("Segoe UI", 8));
    ui->tableViewV->setModel(new QStandardItemModel);
    ui->tableViewV->setFont(QFont("Segoe UI", 8));
    ui->tableViewZ->setModel(new QStandardItemModel);
    ui->tableViewZ->setFont(QFont("Segoe UI", 8));

    ui->progressLabel->hide();
}

MainWindow::~MainWindow()
{
    delete ui;
}

template <typename ValueT, typename HandlerT>
bool MainWindow::checkInput(ValueT &value, QLineEdit *lineEdit, const QString &errorText, HandlerT &handler) {
    bool success = false;
    string inputString;
    try {
        value = handler(lineEdit->text());
        return true;
    }
    catch(const char* error) {
        QMessageBox::critical(this, "Неверное значение!", errorText);
        return false;
    }
}

int64_t intPositiveHandler(const QString &inputString) {
    bool ok;
    int64_t val = inputString.toLongLong(&ok);
    if (ok == false || val <= 0)
        throw "error";
    return val;
}

int64_t intGridHandler(const QString &inputString) {
    bool ok;
    int64_t val = inputString.toLongLong(&ok);
    if (ok == false || val <= 0 || val % 4 != 0)
        throw "error";
    return val;
}

double doublePositiveHandler(const QString &inputString) {
    bool ok;
    double val = inputString.toDouble(&ok);
    if (ok == false || val <= 0.0)
        throw "error";
    return val;
}

double doubleAnyHandler(const QString &inputString) {
    bool ok;
    double val = inputString.toDouble(&ok);
    if (ok == false)
        throw "error";
    return val;
}

double doubleZeroOneHandler(const QString &inputString) {
    bool ok;
    double val = inputString.toDouble(&ok);
    if (ok == false || val < 0.0 || val > 1.0)
        throw "error";
    return val;
}

vector<double> getU(size_t m, size_t n, const vector<size_t> &row_bound) {
    vector<double> uVec((m+1)*(n+1));
    double h = (b-a)/n, k = (d-c)/m;
    size_t N = n+1;

    for (size_t i = 0; i <= m; ++i) {
        double y = c + i*k;
        double x = a + row_bound[2*i]*h;
        for (size_t j = row_bound[2*i]; j < row_bound[2*i+1]; ++j) {
            uVec[N * i + j] = u(x, y);
            x += h;
        }
    }

    return uVec;
}

vector<double> getV0(size_t m, size_t n, const vector<size_t> &row_bound) {
    vector<double> v0Vec((m+1)*(n+1));
    double *v0Ptr = v0Vec.data();
    double h = (b-a)/n, k = (d-c)/m;
    size_t N = n+1;

    int m14 = m/4;
    int m34 = 3*m/4;
    int n14 = n/4;
    int n34 = 3*n/4;

    // горизонтальные границы
    for (size_t i = n/4; i <= 3*n/4; ++i)
        v0Ptr[i] = u(a + i*h, c);
    for (size_t i = 0; i <= n/4; ++i)
        v0Ptr[N*m14 + i] = u(a + i*h, c + m14*k);
    for (size_t i = 3*n/4; i <= n; ++i)
        v0Ptr[N*m14 + i] = u(a + i*h, c + m14*k);
    for (size_t i = 0; i <= n/4; ++i)
        v0Ptr[N*m34 + i] = u(a + i*h, c + m34*k);
    for (size_t i = n/4; i <= n; ++i)
        v0Ptr[N*m + i] = u(a + i*h, d);

    // вертикальные границы
    for (size_t j = m/4; j <= 3*m/4; ++j)
        v0Ptr[N*j] = u(a, c + j*k);
    for (size_t j = 0; j <= m/4; ++j)
        v0Ptr[N*j + n14] = u(a + n14*h, c + j*k);
    for (size_t j = 3*m/4; j <= m; ++j)
        v0Ptr[N*j + n14] = u(a + n14*h, c + j*k);
    for (size_t j = 0; j <= m/4; ++j)
        v0Ptr[N*j + n34] = u(a + n34*h, c + j*k);
    for (size_t j = m/4; j <= m; ++j)
        v0Ptr[N*j + n] = u(b, c + j*k);

    return v0Vec;
}

vector<double> getZ(size_t m, size_t n, const vector<double> &uVec, const vector<double> &vVec, const vector<size_t> &row_bound) {
    vector<double> zVec((m+1)*(n+1));
    size_t N = n+1;

    for (size_t i = 0; i <= m; ++i)
        for (size_t j = row_bound[2*i]; j < row_bound[2*i+1]; ++j)
            zVec[N * i + j] = uVec[N * i + j] - vVec[N * i + j];

    return zVec;
}

void plotSurface(int64_t n, int64_t m, const vector<double> &vec, QSurface3DSeries *series) {
    int64_t N = n+1;
    double h = (b-a)/n, k = (d-c)/m;
    const double *dataPtr = vec.data();
    QSurfaceDataArray *dataArray  = new QSurfaceDataArray;

    for (int64_t i = 0; i <= m; ++i) {
        QSurfaceDataRow *dataRow  = new QSurfaceDataRow;
        const double *currDataPtr = dataPtr + N*i;
        double y = c + i*k;
        double x = a;
        for (int64_t j = 0; j <= n; ++j) {
            *dataRow << QVector3D(x, *(currDataPtr++),  y);
            x += h;
        }
        *dataArray << dataRow;
    }
    series->dataProxy()->resetArray(dataArray);
}

void fillTable(int64_t n, int64_t m, const vector<double> &vec, QStandardItemModel *model) {
    double h = (b-a)/n, k = (d-c)/m;
    double x, y;
    QFont boldFont("Segoe UI", 8);

    boldFont.setBold(true);
    model->setRowCount(m+1);
    model->setColumnCount(n+1);

    x = a;
    for (int j = 0; j <= n; ++j) {
        model->setHeaderData(j, Qt::Horizontal, "x" + QString::number(j) + "=" + QString::number(x, 'f', 3));
        model->horizontalHeaderItem(j)->setFont(boldFont);
        x += h;
    }
    y = c;
    model->setHeaderData(0, Qt::Vertical, "");
    for (int i = 0; i <= m; ++i) {
        model->setHeaderData(i, Qt::Vertical, "y" + QString::number(i) + "=" + QString::number(y, 'f', 3));
        model->verticalHeaderItem(i)->setFont(boldFont);
        y += k;
    }

    int64_t vecIdx = 0;
    for (int64_t i = 0; i <= m; ++i)
        for (int64_t j = 0; j <= n; ++j)
            model->setData(model->index(i, j), QString::number(vec[vecIdx++], 'f', 8));
}


void MainWindow::on_pushButton_clicked()
{
    int64_t m, n, K, nmax;
    double eps;
    if (!checkInput(m, ui->lineEdit_m, "'m' должно быть натуральным числом и кратным 4!", intGridHandler) ||
            !checkInput(n, ui->lineEdit_n, "'n' должно быть натуральным числом и кратным 4!", intGridHandler) ||
            !checkInput(K, ui->lineEdit_k, "'k' должно быть натуральным числом!", intPositiveHandler) ||
            !checkInput(eps, ui->lineEdit_eps, "'eps' должно быть положительным действительным числом!", doublePositiveHandler) ||
            !checkInput(nmax, ui->lineEdit_nmax, "'N_max' должно быть натуральным числом!", intPositiveHandler))
        return;

    ui->progressLabel->show();
    ui->pushButton->setDisabled(true);
    ui->pushButton_2->setDisabled(true);
    ui->checkBox_u->setDisabled(true);
    ui->checkBox_v0->setDisabled(true);
    ui->checkBox_vn->setDisabled(true);
    ui->checkBox_zn->setDisabled(true);

    Worker *worker = new Worker;
    QThread *thread = new QThread;
    worker->moveToThread(thread);
    connect(this, SIGNAL(startMethod(size_t,size_t,size_t,double,size_t)), worker, SLOT(launchMethod(size_t,size_t,size_t,double,size_t)));
    connect(worker, SIGNAL(sendResult(const solution_t&)), this, SLOT(processResult(const solution_t&)));
    connect(worker, SIGNAL(updateProgress(const QString&)), ui->progressLabel, SLOT(setText(const QString&)));
    connect(worker, SIGNAL(finished()), thread, SLOT(quit()));
    connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    thread->start();
    emit startMethod(n, m, K, eps, nmax);
}

void MainWindow::processResult(const solution_t& soln) {
    size_t n = soln.n;
    size_t m = soln.m;
    size_t K = soln.k;
    size_t nmax = soln.iteration_num_max;
    double eps = soln.eps_max;

    size_t N = n+1;
    double h = (b-a)/n, k = (d-c)/m;
    double sinmaxh = sin(M_PI*(n-1)/(2*n));
    double sinmaxk = sin(M_PI*(m-1)/(2*m));
    double Mmax = 4.0 * ((sinmaxh*sinmaxh)/(h*h) + (sinmaxk*sinmaxk)/(k*k)),
           Mmin = 15.0;


    vector<double> uVec = getU(m, n, soln.row_bound);
    vector<double> v0Vec = getV0(m, n, soln.row_bound);
    vector<double> zVec = getZ(m, n, uVec, soln.vVec, soln.row_bound);
    const double *uPtr  = uVec.data(),
                 *vPtr  = soln.vVec.data();

    // ПОДСЧЁТ ДАННЫХ ДЛЯ СПРАВКИ, ЗАПОЛНЕНИЕ СПРАВКИ
    double z_max = -1.0, x_max = -1.0, y_max = -1.0;
    for (size_t j = 1; j < m; ++j) {
        size_t shift = j*N + soln.row_bound[2*j];
        const double *uc = uPtr + shift;
        const double *vc = vPtr + shift;
        for (size_t i = soln.row_bound[2*j]; i < soln.row_bound[2*j + 1]; ++i) {
            if (abs(*uc - *vc) > z_max) {
                z_max = abs(*uc - *vc);
                x_max = a + h*i;
                y_max = c + k*j;
            }
            ++uc;
            ++vc;
        }
    }


    solutionInfo->hide();

    solutionLabel->setText(QString("Для решения тестовой задачи использована сетка-основа\n\
с числом разбиений по 'x' n = %1 и числом разбиений по 'y' m = %2.\n\n\
Метод - с чебышёвским набором параметров;\n\
параметры: k = %3, M_max = %4, M_min = %5\n\
(расчитаны на основе оценки собственных чисел).\n\
Критерии остановки:\n\
по точности - 'ε_точн = %6',\n\
по числу итераций - 'N_max = %7'.\n\n\
На решение СЛАУ затрачено 'N = %8' итераций,\n\
достигнута точность итерационного метода 'ε_N = %9'.\n\n\
Невязка СЛАУ на решении имеет максимум-норму ||R_N||∞ = %10\n\
и норму ||R_N||₂ = %15.\n\n\
Тестовая задача должна быть решена с погрешностью не более ε₀ = 0.5e-6,\n\
а решена с погрешностью ε₁ = %11.\n\
Максимальное отклонение точного и численного решений\n\
наблюдается в узле (x,y) = (%12,%13).\n\n\
В качестве начального приближения взят нулевой вектор.\n\n\
Невязка СЛАУ на начальном приближении имеет евклидову норму ||R₀||₂ = %14.")
            .arg(n).arg(m).arg(K).arg(Mmax).arg(Mmin).arg(eps).arg(nmax).arg(soln.iteration_num).arg(soln.eps)
            .arg(soln.r_max).arg(z_max).arg(x_max).arg(y_max).arg(QString::number(soln.r0_max, 'g', 10)).arg(soln.r_2));

    solutionInfo->show();
    solutionInfo->showNormal();

    // ПОСТРОЕНИЕ ГРАФИКОВ
    plotSurface(n, m, uVec, uSeries);
    plotSurface(n, m, v0Vec, v0Series);
    plotSurface(n, m, soln.vVec, vnSeries);
    plotSurface(n, m, zVec, znSeries);

    // ЗАПОЛНЕНИЕ ТАБЛИЦ
    QStandardItemModel *uModel = (QStandardItemModel*)ui->tableViewU->model();
    QStandardItemModel *vModel = (QStandardItemModel*)ui->tableViewV->model();
    QStandardItemModel *zModel = (QStandardItemModel*)ui->tableViewZ->model();
    fillTable(n, m, uVec, uModel);
    fillTable(n, m, soln.vVec, vModel);
    fillTable(n, m, zVec, zModel);

    ui->progressLabel->hide();
    ui->pushButton->setEnabled(true);
    ui->pushButton_2->setEnabled(true);
    ui->checkBox_u->setEnabled(true);
    ui->checkBox_v0->setEnabled(true);
    ui->checkBox_vn->setEnabled(true);
    ui->checkBox_zn->setEnabled(true);
}

void MainWindow::on_checkBox_u_clicked()
{
    if (ui->checkBox_u->isChecked())
        graph->addSeries(uSeries);
    else
        graph->removeSeries(uSeries);
}


void MainWindow::on_checkBox_v0_clicked()
{
    if (ui->checkBox_v0->isChecked())
        graph->addSeries(v0Series);
    else
        graph->removeSeries(v0Series);
}


void MainWindow::on_checkBox_vn_clicked()
{
    if (ui->checkBox_vn->isChecked())
        graph->addSeries(vnSeries);
    else
        graph->removeSeries(vnSeries);
}


void MainWindow::on_checkBox_zn_clicked()
{
    if (ui->checkBox_zn->isChecked())
        graph->addSeries(znSeries);
    else
        graph->removeSeries(znSeries);
}


void MainWindow::on_pushButton_2_clicked()
{
    if (!solutionInfo->isHidden())
        solutionInfo->hide();
    solutionInfo->show();
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
    int key = event->key();

    switch (key) {
        case Qt::Key_Enter:
        case Qt::Key_Return:
            ui->pushButton->click();
            break;
        case Qt::Key_Up:
            if (ui->lineEdit_nmax->hasFocus()) {
                ui->lineEdit_nmax->clearFocus();
                ui->lineEdit_eps->setFocus();
            } else if (ui->lineEdit_eps->hasFocus()) {
                ui->lineEdit_eps->clearFocus();
                ui->lineEdit_k->setFocus();
            } else if (ui->lineEdit_k->hasFocus()) {
                ui->lineEdit_k->clearFocus();
                ui->lineEdit_n->setFocus();
            } else if (ui->lineEdit_n->hasFocus()) {
                ui->lineEdit_n->clearFocus();
                ui->lineEdit_m->setFocus();
            } else if (ui->lineEdit_m->hasFocus()) {
                ui->lineEdit_m->clearFocus();
                ui->lineEdit_nmax->setFocus();
            } else
                ui->lineEdit_nmax->setFocus();
            break;
        case Qt::Key_Down:
            if (ui->lineEdit_m->hasFocus()) {
                ui->lineEdit_m->clearFocus();
                ui->lineEdit_n->setFocus();
            } else if (ui->lineEdit_n->hasFocus()) {
                ui->lineEdit_n->clearFocus();
                ui->lineEdit_k->setFocus();
            } else if (ui->lineEdit_k->hasFocus()) {
                ui->lineEdit_k->clearFocus();
                ui->lineEdit_eps->setFocus();
            } else if (ui->lineEdit_eps->hasFocus()) {
                ui->lineEdit_eps->clearFocus();
                ui->lineEdit_nmax->setFocus();
            } else if (ui->lineEdit_nmax->hasFocus()) {
                ui->lineEdit_nmax->clearFocus();
                ui->lineEdit_m->setFocus();
            } else
                ui->lineEdit_m->setFocus();
            break;
    }
}
