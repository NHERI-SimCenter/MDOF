#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <math.h>

class MyGlWidget;
class Vector;
class EarthquakeRecord;

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
    float setSelectionBoundary(float y1, float y2);

private slots:
    void on_inFloors_editingFinished();
    void on_inWeight_editingFinished();
    void on_inHeight_editingFinished();

    void on_stopButton_clicked();
    void on_runButton_clicked();
 //   void on_slider_sliderMoved(int position);

    void on_slider_valueChanged(int value);
    void on_slider_sliderPressed();
    void on_slider_sliderReleased();

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
    int numSteps;
    double *gMotion;
    Vector *elCentroData;

    bool needAnalysis;
    double **dispResponses;
    double maxDisp;
    int currentStep;
    bool stopRun;

    bool movingSlider;
    int fMinSelected, fMaxSelected;
    int sMinSelected, sMaxSelected;
};

#endif // MAINWINDOW_H
