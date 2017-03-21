#ifndef MAINWINDOW_H
#define MAINWINDOW_H

/* *****************************************************************************
Copyright (c) 2016-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: fmckenna


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
    Vector *eqData;

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
