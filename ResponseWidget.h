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

class ResponseWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ResponseWidget(MainWindow *main, int mainWindowItem, QString &label, QString &xAxis, QString &yAxis, QWidget *parent = 0);
    ~ResponseWidget();

    int getItem();
    void setItem(int);
    void setData(QVector<double> &data, QVector<double> &time, int numSteps, double dt);
    void setData(QVector<double> &data, QVector<double> &x, int numSteps);

signals:

public slots:
    void itemEditChanged(void);

private:
    QCustomPlot *thePlot;
    QLineEdit *theItemEdit;
    int theItem; // floor or story #
    int mainWindowItem; // tag used in call back

    MainWindow *main;

    QCPGraph *graph;
    QCPItemTracer *groupTracer;
    QCPCurve *curve;

};

#endif // NODERESPONSEWIDGET_H
