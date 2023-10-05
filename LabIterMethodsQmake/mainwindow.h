#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <Q3DSurface>
#include <QLineEdit>
#include <QLabel>
#include "worker.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    Q3DSurface *graph;
    QSurface3DSeries *uSeries;
    QSurface3DSeries *v0Series;
    QSurface3DSeries *vnSeries;
    QSurface3DSeries *znSeries;
    ~MainWindow();
public slots:
    void processResult(const solution_t&);

signals:
    void startMethod(size_t, size_t, size_t, double, size_t);

private slots:
    void on_pushButton_clicked();

    void on_checkBox_u_clicked();

    void on_checkBox_v0_clicked();

    void on_checkBox_vn_clicked();

    void on_checkBox_zn_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;
    QWidget *plot;
    QWidget *solutionInfo;
    QLabel *solutionLabel;

    bool uActive = false,
         v0Active = false,
         vnActive = false,
         znActive = false;
    template <typename ValueT, typename HandlerT>
    bool checkInput(ValueT &value, QLineEdit *lineEdit, const QString &errorText, HandlerT &handler);
    void keyPressEvent(QKeyEvent *event);
};

#endif // MAINWINDOW_H
