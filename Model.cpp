#include "Model.h"
#include <QDebug>

Model::Model(QObject *parent) :
    QObject(parent)
{
}

  void
  Model::setNumFloors(int num){
      numFloors = num;
      if (mass != 0 || numFloors != num) {
          delete [] mass;
      }
      mass = new double[num];
      qDebug() << numFloors;
  }

  void Model::setMasses(double *newMasses, int floor)
  {
      if (floor != 0)
           mass[floor-1] = newMasses[0];
      else {
           for (int i=0; i<numFloors; i++)
               mass[i]=newMasses[i];
        }

  }
