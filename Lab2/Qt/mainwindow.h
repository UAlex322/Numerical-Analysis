#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLabel>


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QLabel  *windowLabel;

private slots:
    void on_tButtonSolve_pressed();

    void on_mButtonSolve_pressed();

    void on_tButtonSolnInfo_clicked();

private:
    QWidget *testSolutionInfo,
            *window;
    Ui::MainWindow *ui;
    bool checkInput(const QString &, int &);
    void InfoWidget(QWidget *parent, QLabel *label);
};
#endif // MAINWINDOW_H
