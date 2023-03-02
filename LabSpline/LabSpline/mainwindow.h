#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <vector>
#include <functional>
#include "spline.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_buttonSolve_clicked();

    void on_comboBoxFunction_activated(int index);

    void on_pushButton_2_clicked();

    void on_checkBoxF_clicked();

    void on_checkBoxS_clicked();

    void on_checkBoxF1_clicked();

    void on_checkBoxS1_clicked();

    void on_checkBoxF2_clicked();

    void on_checkBoxS2_clicked();

    void on_checkBoxD_clicked();

    void on_checkBoxD1_clicked();

    void on_checkBoxD2_clicked();

    void on_buttonSolnInfo_clicked();

private:
    Ui::MainWindow *ui;
    QWidget *info,
            *solutionInfo;
    QLabel *infoLabel,
           *solutionInfoLabel;

    double t;
    std::function<double(double)> func[3];
    std::function<double(double)> funcD[3];
    std::function<double(double)> funcD2[3];
    std::vector<double> x;
};
#endif // MAINWINDOW_H
