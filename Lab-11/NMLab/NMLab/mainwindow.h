#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QStandardItemModel>
#include <QTabWidget>
#include <QTableView>
#include <QRadioButton>
#include <QGraphicsScene>
#include <QLabel>
#include <qcustomplot.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_progInfoButton_clicked();

    void on_main2SolveButton_clicked();

    void on_main2ShowSolnInfoButton_clicked();

    void on_tableInfoButton_pressed();

    void on_clearButton_pressed();

private:
    Ui::MainWindow *ui;

    QWidget *main2SolutionInfo,
            *progInfo,
            *tableInfo;
    QLabel *main2SolutionInfoLabel,
           *progInfoLabel,
           *tableInfoLabel;
    std::vector<QCPCurve*> phase;

    bool main2IsSolvedOnce = false;
    size_t graphCount = 0;

    void initializeOutputInfoWidget(QWidget *parent, QLabel *label);
    void initializePlotWidget(QCustomPlot *plotWidget, QWidget *backgroundWidget);
    void initializeBackground(QWidget *widget, QWidget *backgroundWidget);
    void initializeMainTaskTable(QTableView *table, size_t rows);
    void setGraphLabels(double L, double x0, double u0, double ud0);
    bool checkInput(const QString &, const QString &, const QString &, const QString &, const QString &, const QString &, const QString &,
                    double &, double &, double &, double &, double &, int &, double &);
    bool checkMain2Input(const QString &, const QString &, double &, double &);
};

#endif // MAINWINDOW_H
