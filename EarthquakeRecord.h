#ifndef EARTHQUAKERECORD_H
#define EARTHQUAKERECORD_H

#include <QString>
class Vector;

class EarthquakeRecord
{
public:
    EarthquakeRecord(QString fileName);
    EarthquakeRecord(QString name, int numSteps, double dt, Vector *data);
    ~EarthquakeRecord();

    QString name;
    int numSteps;
    double dt;
    Vector *data;
};

#endif // EARTHQUAKERECORD_H
