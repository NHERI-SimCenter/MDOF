#ifndef MODEL_H
#define MODEL_H

#include <QObject>

class Model : public QObject
{
    Q_OBJECT

public:
    explicit Model(QObject *parent = 0);

    setNumFloors(int numFloors);

    setMasses(double *newMasses, int floor = 0);
    setStoryData(double *stiff,double *fy, double *b, double *h, int story=0);
    setMotion(int numData, double dT, double *data);

    analyze();

    draw(int step);
    getMaxResponse(char *response);

signals:

public slots:

    private:
    int numFloors;
    
    double *mass;
    double *stiffness;
    double *fy;
    double *height
};

#endif // MODEL_H
