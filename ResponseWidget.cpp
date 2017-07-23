#include "ResponseWidget.h"
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



ResponseWidget::ResponseWidget(MainWindow *mainWindow, int mainItem, QWidget *parent)
    : QWidget(parent), theItem(0), mainWindowItem(mainItem), main(mainWindow)
{
    // create a main layout
    QVBoxLayout *mainLayout = new QVBoxLayout();

    theItemEdit = createTextEntry(tr("Item"), mainLayout);
    theItemEdit->setValidator(new QIntValidator);
    connect(theItemEdit, SIGNAL(editingFinished()), this, SLOT(itemEditChanged()));

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
    theItemEdit->setText(QString::number(newItem));
}

void
ResponseWidget::itemEditChanged() {
    int oldItem = theItem;
    QString textItems =  theItemEdit->text();
    theItem = textItems.toInt();
    if (theItem <= 0 || theItem > main->getNumFloors()) {
        theItem=oldItem;
        theItemEdit->setText(QString::number(theItem));
    } else {
        main->setResponse(theItem, mainWindowItem);
    }
}

void
ResponseWidget::setData(QVector<double> &data, QVector<double> time, int numSteps, double dt) {
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

void
ResponseWidget::setData(QVector<double> &data, QVector<double> x, int numSteps) {
    thePlot->clearGraphs();
    graph = thePlot->addGraph();

    thePlot->graph(0)->setData(x, data, true);

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
        qDebug() << i*.02 << " " << value << " " << xValue;
    }

   //
   thePlot->yAxis->setRange(minValue, maxValue);
   thePlot->xAxis->setRange(xMinValue, xMaxValue);

    //thePlot->axisRect()->setAutoMargins(QCP::msNone);
   // thePlot->axisRect()->setMargins(QMargins(0,0,0,0));
    thePlot->replot();

}
