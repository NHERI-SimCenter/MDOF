#include "EarthquakeRecord.h"
#include <Vector.h>

EarthquakeRecord::EarthquakeRecord(QString fileName)
{

}

EarthquakeRecord::EarthquakeRecord(QString theName, int numberSteps, double theDt, Vector *theData)
    :name(theName), numSteps(numberSteps), dt(theDt), data(0)
{
     data = new Vector(*theData);
}

EarthquakeRecord::~EarthquakeRecord()
{
    if (data != 0)
        delete [] data;
}

