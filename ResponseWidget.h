#ifndef RESPONSEWIDGET_H
#define RESPONSEWIDGET_H

#include <QWidget>
#include <QVector>

class QCustomPlot;
class QLineEdit;
class QCPGraph;
class QCPCurve;
class QCPItemTracer;
class MainWindow;
class QSpinBox;


//
// a widget to plot the response,
//  NOTE to cut down on memory currently the widget does not keep data, requires a call back to main window
//  THIS is proving problematic to say the least! .. will redo
//

class ResponseWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ResponseWidget(MainWindow *main, int mainWindowItem, QString &label, QString &xAxis, QString &yAxis, QWidget *parent = 0);
    ~ResponseWidget();

    int getItem();
    void setItem(int);
    void setData(QVector<double> &data, QVector<double> &time, int numSteps, double dt, int index);
    void setData(QVector<double> &data, QVector<double> &x, int numSteps, int index);

signals:

public slots:
    void itemEditChanged(int);

private:
    QCustomPlot *thePlot;
    QSpinBox *theSpinBox;

    int theItem; // floor or story #

    // 3 variables neeeded for call back .. need to go away!
    int mainWindowItem; 
    MainWindow *main;   
    bool dataSET;

    QCPGraph *graph;
    QCPItemTracer *groupTracer;
    QCPCurve *curve;
};

#endif // NODERESPONSEWIDGET_H
