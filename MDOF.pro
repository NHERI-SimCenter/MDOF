#-------------------------------------------------
#
# Project created by QtCreator 2017-03-11T14:00:13
#
#-------------------------------------------------

QT       += core gui opengl printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = MDOF
TEMPLATE = app

include(OPS_includes.pro)


SOURCES += main.cpp\
        MainWindow.cpp \
    MyGlWidget.cpp \
    qcustomplot.cpp \
    EarthquakeRecord.cpp

HEADERS  += MainWindow.h \
    MyGlWidget.h \
    qcustomplot.h \
    EarthquakeRecord.h

FORMS    += MainWindow.ui
