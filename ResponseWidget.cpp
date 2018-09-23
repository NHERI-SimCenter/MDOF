#include "ResponseWidget.h"
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <MainWindow.h>
#include <QSpinBox>

#include <qcustomplot.h>


ResponseWidget::ResponseWidget(MainWindow *mainWindow,
                               int mainItem,
                               QString &label,
                               QString &xLabel,
                               QString &yLabel,
                               QWidget *parent)
    : QWidget(parent), theItem(0), mainWindowItem(mainItem), main(mainWindow)
{
    // create a main layout
    QVBoxLayout *mainLayout = new QVBoxLayout();


    //
    // spin box
    //

    dataSET = false;
    theSpinBox = new QSpinBox();
    theSpinBox->setMinimum(1);
    theSpinBox->setMaximum(99);
    theSpinBox->setMinimumWidth(100);
    QHBoxLayout *theSpinLayout = new QHBoxLayout();
    QLabel *spinLabel = new QLabel(label);
    theSpinLayout->addWidget(spinLabel);
    theSpinLayout->addWidget(theSpinBox);
    theSpinLayout->addStretch();
    connect(theSpinBox,SIGNAL(valueChanged(int)), this, SLOT(itemEditChanged(int)));
    mainLayout->addLayout(theSpinLayout);

    //
    // graphic window
    //

    thePlot=new QCustomPlot();
    thePlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);

    QRect rec = QApplication::desktop()->screenGeometry();

    int height = 0.2*rec.height();
    int width = 0.5*rec.width();

    thePlot->setMinimumWidth(width);
    thePlot->setMinimumHeight(height);
    mainLayout->addWidget(thePlot);

    thePlot->xAxis->setLabel(xLabel);
    thePlot->yAxis->setLabel(yLabel);

    this->setLayout(mainLayout);
}

ResponseWidget::~ResponseWidget()
{

}

int
ResponseWidget::getItem() {
    return theItem;
}

void
ResponseWidget::setItem(int newItem) {
    theItem = newItem;
    theSpinBox->setValue(newItem);
}

void
ResponseWidget::itemEditChanged(int theItem) {
    if (dataSET == true)
        main->setResponse(theItem, mainWindowItem);
}


void
ResponseWidget::setData(QVector<double> &data, QVector<double> &time, int numSteps, double dt, int index = -1) {

    dataSET = false;
    if (index != -1) {
        theSpinBox->setSingleStep(1);
        theSpinBox->setValue(index);

    }
    dataSET = true;

    if (time.size() != data.size()) {
        qDebug() << "ResponseWidget - setData vectors of differing sizes";
        return;
    }

    if (time.size() != numSteps) {
        qDebug() << "ResponseWidget - setData vector and step size do not agree";
        return;
    }


    thePlot->clearGraphs();
    graph = thePlot->addGraph();

    thePlot->graph(0)->setData(time, data);

    double minValue = 0;
    double maxValue = 0;
    for (int i=0; i<numSteps; i++) {
        double value = data.at(i);
        if (value < minValue)
            minValue = value;
        if (value > maxValue)
            maxValue = value;
    }

    thePlot->xAxis->setRange(0, numSteps*dt);
    thePlot->yAxis->setRange(minValue, maxValue);
    //thePlot->axisRect()->setAutoMargins(QCP::msNone);
    thePlot->axisRect()->setMargins(QMargins(0,0,0,0));
    thePlot->replot();

}

void
ResponseWidget::setData(QVector<double> &data, QVector<double> &x, int numSteps, int index) {

    dataSET = false;
    if (index != -1) {
        theSpinBox->setRange(1, main->getNumFloors());
        theSpinBox->setValue(index);
    }
    dataSET = true;

    thePlot->clearGraphs();
    thePlot->clearPlottables();
    curve = new QCPCurve(thePlot->xAxis, thePlot->yAxis);

    curve->setData(x,data);

    //thePlot->graph(0)->setData(x, data, true);

    double minValue = 0;
    double maxValue = 0;
    double xMinValue = 0;
    double xMaxValue = 0;
    for (int i=0; i<numSteps; i++) {
        double value = data.at(i);
        double xValue = x.at(i);
        if (value < minValue)
            minValue = value;
        if (value > maxValue)
            maxValue = value;
        if (xValue < xMinValue)
            xMinValue = xValue;
        if (xValue > xMaxValue)
            xMaxValue = xValue;
    }

   thePlot->yAxis->setRange(minValue, maxValue);
   thePlot->xAxis->setRange(xMinValue, xMaxValue);

    //thePlot->axisRect()->setAutoMargins(QCP::msNone);
   // thePlot->axisRect()->setMargins(QMargins(0,0,0,0));
    thePlot->replot();
}
