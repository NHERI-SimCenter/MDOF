#include "NodeResponseWidget.h"
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <MainWindow.h>

#include <qcustomplot.h>

extern QLineEdit *
createTextEntry(QString text,
                QVBoxLayout *theLayout,
                int minL=100,
                int maxL=100);



NodeResponseWidget::NodeResponseWidget(MainWindow *mainWindow, QWidget *parent)
    : QWidget(parent), theFloor(0), main(mainWindow)
{
    // create a main layout
    QVBoxLayout *mainLayout = new QVBoxLayout();

    theFloorEdit = createTextEntry(tr("Floor"), mainLayout);
    theFloorEdit->setValidator(new QIntValidator);
    connect(theFloorEdit, SIGNAL(editingFinished()), this, SLOT(floorEditChanged()));

    thePlot=new QCustomPlot();
    thePlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);

    QRect rec = QApplication::desktop()->screenGeometry();

    int height = 0.2*rec.height();
    int width = 0.5*rec.width();

    thePlot->setMinimumWidth(width);
    thePlot->setMinimumHeight(height);
    mainLayout->addWidget(thePlot);

    this->setLayout(mainLayout);
}

NodeResponseWidget::~NodeResponseWidget()
{

}

int
NodeResponseWidget::getFloor() {
    return theFloor;
}

void
NodeResponseWidget::setFloor(int newFloor) {
    theFloor = newFloor;
    theFloorEdit->setText(QString::number(newFloor));
}

void
NodeResponseWidget::floorEditChanged() {
    int oldFloor = theFloor;
    QString textFloors =  theFloorEdit->text();
    theFloor = textFloors.toInt();
    if (theFloor < 0 || theFloor > main->getNumFloors()) {
        theFloor=oldFloor;
        theFloorEdit->setText(QString::number(theFloor));
    } else {
        main->setFloorResponse(theFloor);
    }
}

void
NodeResponseWidget::setData(QVector<double> &data, QVector<double> time, int numSteps, double dt) {
    thePlot->clearGraphs();
    graph = thePlot->addGraph();
    thePlot->graph(0)->setData(time, data);
    thePlot->xAxis->setRange(0, numSteps*dt);
    double minValue = 0;
    double maxValue = 0;
    for (int i=0; i<numSteps; i++) {
        double value = data.at(i);
        if (value < minValue)
            minValue = value;
        if (value > maxValue)
            maxValue = value;
    }

    thePlot->yAxis->setRange(minValue, maxValue);
    //thePlot->axisRect()->setAutoMargins(QCP::msNone);
    thePlot->axisRect()->setMargins(QMargins(0,0,0,0));
    thePlot->replot();

}
