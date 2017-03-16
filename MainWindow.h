#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <math.h>

class MyGlWidget;
class Vector;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void draw(MyGlWidget *);
    void doAnalysis();
    float getHeight() {return buildingH;};
    float getMaxDisp(){return maxDisp;};

private slots:
    void on_inFloors_editingFinished();
    void on_inWeight_editingFinished();
    void on_inHeight_editingFinished();

    void on_run_pressed();

private:
    void updatePeriod();
    void setBasicModel(int numFloors, double period);
    void setBasicModel(int numFloors, double buildingWeight, double buildingK);
    

private:
    Ui::MainWindow *ui;

    int numFloors;
    double period;
    double buildingW;
    double buildingH;
    double storyK;

    double *masses;
    double *k;
    double *fy;
    double *b;
    double *heights;
    double dampingRatio;

    double dt;
    double *gMotion;

    bool needAnalysis;
    Vector *elCentroData;
    double **dispResponses;
    double maxDisp;
    int currentTime;
};

#endif // MAINWINDOW_H
