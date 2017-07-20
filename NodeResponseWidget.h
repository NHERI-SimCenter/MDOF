#ifndef NODERESPONSEWIDGET_H
#define NODERESPONSEWIDGET_H

#include <QWidget>
#include <QVector>

class QCustomPlot;
class QLineEdit;
class QCPGraph;
class QCPItemTracer;
class MainWindow;

class NodeResponseWidget : public QWidget
{
    Q_OBJECT
public:
    explicit NodeResponseWidget(MainWindow *main, QWidget *parent = 0);
    ~NodeResponseWidget();

    int getFloor();
    void setFloor(int);
    void setData(QVector<double> &data, QVector<double> time, int numSteps, double dt);

signals:

public slots:
    void floorEditChanged(void);

private:
    QCustomPlot *thePlot;
    QLineEdit *theFloorEdit;
    int theFloor;

    MainWindow *main;

    QCPGraph *graph;
    QCPItemTracer *groupTracer;
};

#endif // NODERESPONSEWIDGET_H
