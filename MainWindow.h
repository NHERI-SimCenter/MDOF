#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//#include "PropertiesWidget.h"
#include <QMainWindow>
#include <math.h>
#include <map>
#include <QString>

class MyGlWidget;
class Vector;
class EarthquakeRecord;
class QCPGraph;
class QCPItemTracer;
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

    friend class PropertiesWidget;


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

    void on_inFloorWeight_editingFinished();
 //   void on_inFloorHeight_editingFinished();

    void on_inStoryHeight_editingFinished();
    void on_inStoryK_editingFinished();
    void on_inStoryFy_editingFinished();
    void on_inStoryB_editingFinished();

    void on_tableWidget_cellChanged(int row, int column);
    void on_tableWidget_cellClicked(int row, int column);

    void on_inGravity_editingFinished();

    void on_inMotionSelection_currentTextChanged(const QString &arg1);

private:
    void updatePeriod();
    void setBasicModel(int numFloors, double period);
    void setBasicModel(int numFloors, double buildingWeight, double buildingK);
    void reset(void);
    
private:
    Ui::MainWindow *ui;

    int numFloors;
    double period;
    double buildingW;
    double buildingH;
    double storyK;

    double *weights;
    double *k;
    double *fy;
    double *b;
    double *floorHeights;
    double *storyHeights;
    double dampingRatio;
    double g;

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

    bool updatingPropertiesTable;

    QVector<double> time;
    QVector<double> values;
     QCPGraph *graph;
     QCPItemTracer *groupTracer;

     std::map <QString, EarthquakeRecord *> records;
};

#endif // MAINWINDOW_H
