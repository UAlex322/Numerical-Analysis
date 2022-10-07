#-------------------------------------------------
#
# Project created by QtCreator 2022-09-26T00:59:32
#
#-------------------------------------------------

QT       += core gui printsupport
CONFIG += c++17
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = NMLab
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
        qcustomplot.cpp

HEADERS  += mainwindow.h \
    Methods.h \
    qcustomplot.h

FORMS    += mainwindow.ui
