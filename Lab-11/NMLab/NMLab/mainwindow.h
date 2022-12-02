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

private:
    Ui::MainWindow *ui;

    QWidget *testSolutionInfo,
            *main1SolutionInfo,
            *main2SolutionInfo,
            *progInfo,
            *tableInfo;
    QLabel *testSolutionInfoLabel,
           *main1SolutionInfoLabel,
           *main2SolutionInfoLabel,
           *progInfoLabel,
           *tableInfoLabel;
    QCPCurve *phase;

    bool testIsSolvedOnce = false;
    bool main1IsSolvedOnce = false;
    bool main2IsSolvedOnce = false;

    void initializeOutputInfoWidget(QWidget *parent, QLabel *label);
    void initializePlotWidget(QCustomPlot *plotWidget, QWidget *backgroundWidget);
    void initializeBackground(QWidget *widget, QWidget *backgroundWidget);

    //void initializeTestTaskTable(QTableView *table, size_t rows);
    void initializeMainTaskTable(QTableView *table, size_t rows);
    bool checkInput(const QString &, const QString &, const QString &, const QString &, const QString &, const QString &, const QString &,
                    double &, double &, double &, double &, double &, int &, double &);
    bool checkMain2Input(const QString &, const QString &, double &, double &);
};

#endif // MAINWINDOW_H
