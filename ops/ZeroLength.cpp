/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.27 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLength.cpp,v $

// Written: GLF
// Created: 12/99
// Revision: A
//
// Description: This file contains the implementation for the ZeroLength class.
//
// What: "@(#) ZeroLength.C, revA"

#include "ZeroLength.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>
#include <elementAPI.h>
#include <vector>

// initialise the class wide variables
Matrix ZeroLength::ZeroLengthM2(2,2);
Vector ZeroLength::ZeroLengthV2(2);

void* OPS_ZeroLength()
{
  return 0;
}


//  Construct element with one unidirectional material (numMaterials1d=1)
ZeroLength::ZeroLength(int tag,
		       int Nd1, int Nd2, 
		       UniaxialMaterial &theMat,
		       double pdivl)
 :Element(tag,ELE_TAG_ZeroLength),     
  connectedExternalNodes(2),
  theMatrix(0), theVector(0)
{
  // allocate memory for numMaterials1d uniaxial material models
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
    opserr << "FATAL ZeroLength::ZeroLength - failed to get a copy of material " << theMat.getTag() << endln;
    exit(-1);
  }
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  //opserr << tag << " " << connectedExternalNodes;
  PdivL = pdivl;
}

ZeroLength::ZeroLength(void)
  :Element(0,ELE_TAG_ZeroLength),     
   connectedExternalNodes(2),
   theMatrix(0), theVector(0),
   theMaterial(0)
{
    // ensure the connectedExternalNode ID is of correct size 
    if (connectedExternalNodes.Size() != 2)
      opserr << "FATAL ZeroLength::ZeroLength - failed to create an ID of correct size\n";

    PdivL = 0.;
}


//  Destructor:
//  delete must be invoked on any objects created by the object
//  and on the matertial object.
ZeroLength::~ZeroLength()
{
  delete theMaterial;
}


int
ZeroLength::getNumExternalNodes(void) const
{
    return 2;
}


const ID &
ZeroLength::getExternalNodes(void) 
{
    return connectedExternalNodes;
}



Node **
ZeroLength::getNodePtrs(void) 
{
  return theNodes;
}

int
ZeroLength::getNumDOF(void) 
{
  return 2;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLength element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
ZeroLength::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	return;
    }

    // set default values for error conditions
    numDOF = 2;
    theMatrix = &ZeroLengthM2;
    theVector = &ZeroLengthV2;
    
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if ( theNodes[0] == 0 || theNodes[1] == 0 ) {
      if (theNodes[0] == 0) 
        opserr << "WARNING ZeroLength::setDomain() - Nd1: " << Nd1 << " does not exist in ";
      else
        opserr << "WARNING ZeroLength::setDomain() - Nd2: " << Nd2 << " does not exist in ";

      opserr << "model for ZeroLength ele: " << this->getTag() << endln;

      return;
    }

    this->DomainComponent::setDomain(theDomain);
}   	 



int
ZeroLength::commitState()
{
    int code=0;

    // call element commitState to do any base class stuff
    if ((code = this->Element::commitState()) != 0) {
      opserr << "ZeroLength::commitState () - failed in base class";
    }    

    return theMaterial->commitState();
}

int
ZeroLength::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();
}


int
ZeroLength::revertToStart()
{   
  return theMaterial->revertToStart();
}


int
ZeroLength::update(void)
{
    double strain;
    double strainRate;

    // get trial displacements and take difference
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    Vector  diff  = disp2-disp1;
    const Vector& vel1  = theNodes[0]->getTrialVel();
    const Vector& vel2  = theNodes[1]->getTrialVel();
    Vector  diffv = vel2-vel1;
    
    strain = diff(0);
    strainRate = diffv(0);

    return theMaterial->setTrialStrain(strain,strainRate);
}

const Matrix &
ZeroLength::getTangentStiff(void)
{


    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    
    // get tangent for material
    double E = theMaterial->getTangent();
    double Kg = PdivL;
    double K = E + Kg;

    stiff(0,0) =  K; stiff(0,1) = -K;
    stiff(1,0) = -K; stiff(1,1) =  K;
    // opserr << this->getTag() << " " << K << endln;
    return stiff;
}


const Matrix &
ZeroLength::getInitialStiff(void)
{
    // stiff is a reference to the matrix holding the stiffness matrix
    Matrix& stiff = *theMatrix;
    
    // zero stiffness matrix
    stiff.Zero();
    
    double E = theMaterial->getInitialTangent();
    double Kg = PdivL;
    double K = E + Kg;

    stiff(0,0) =  K; stiff(0,1) = -K;
    stiff(1,0) = -K; stiff(1,1) =  K;

    return stiff;
}
    

const Matrix &
ZeroLength::getDamp(void)
{
    // damp is a reference to the matrix holding the damping matrix
    Matrix& damp = *theMatrix;

    // zero damping matrix
    damp.Zero();

    // get Rayleigh damping matrix 
    int useRayleighDamping = 0;

    if (useRayleighDamping == 1) {

        damp = this->Element::getDamp();
    }
    return damp;
}


const Matrix &
ZeroLength::getMass(void)
{
  // no mass 
  theMatrix->Zero();    
  return *theMatrix; 
}


void 
ZeroLength::zeroLoad(void)
{
  // does nothing now
}

int 
ZeroLength::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ZeroLength::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}

int 
ZeroLength::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}


const Vector &
ZeroLength::getResistingForce()
{
  // zero the residual
  theVector->Zero();
  double force = theMaterial->getStress() + PdivL*theMaterial->getStrain();
  (*theVector)(0) = -force;
  (*theVector)(1) = force;
  return *theVector;
}


const Vector &
ZeroLength::getResistingForceIncInertia()
{	
  // this already includes damping forces from materials
  this->getResistingForce();
  
  // add the damping forces from rayleigh damping
  int useRayleighDamping = 0;
  if (useRayleighDamping == 1) {
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
      *theVector += this->getRayleighDampingForces();
    }  
  } 

  return *theVector;
}


int
ZeroLength::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
ZeroLength::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


int
ZeroLength::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return -1;
}


void
ZeroLength::Print(OPS_Stream &s, int flag)
{
  s << "ZeroLength: " << this->getTag() << "\n nodes: " << this->connectedExternalNodes(0);
  s << " " << connectedExternalNodes(1) << " ";
  s << this->theMaterial->getTangent() << "\n";
}

Response*
ZeroLength::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

      theResponse = new ElementResponse(this, 1, Vector(numDOF));


    } else if ((strcmp(argv[0],"defoANDforce") == 0) ||
        (strcmp(argv[0],"deformationANDforces") == 0) ||
        (strcmp(argv[0],"deformationsANDforces") == 0)) {
      
      theResponse = new ElementResponse(this, 4, Vector(2));
      
    }

    output.endTag();

    return theResponse;
}

int 
ZeroLength::getResponse(int responseID, Information &eleInformation)
{
    const Vector& disp1 = theNodes[0]->getTrialDisp();
    const Vector& disp2 = theNodes[1]->getTrialDisp();
    const Vector  diff  = disp2-disp1;

    switch (responseID) {
    case -1:
        return -1;

    case 1:
        return eleInformation.setVector(this->getResistingForce());


    case 4:
        if (eleInformation.theVector != 0) {
	  (*(eleInformation.theVector))(0) = theMaterial->getStrain();
          (*(eleInformation.theVector))(1) = theMaterial->getStress();
        }
        return 0;      

    default:
        return -1;
    }
}

int
ZeroLength::setParameter(const char **argv, int argc, Parameter &param)
{
  int result = -1;  

  if (argc < 1)
    return -1;

  if (strcmp(argv[0], "material") == 0) {
      if (argc > 2) {
	return theMaterial->setParameter(&argv[1], argc-1, param);
      } else {
	return -1;
      }
  }

  int res = theMaterial->setParameter(argv, argc, param);
  if (res != -1) {
    result = res;
  }

  return result;
}

double
ZeroLength::getForce() {

    double force = theMaterial->getStress() + PdivL*theMaterial->getStrain();

    return force;
}

double
ZeroLength::getDrift(){
    return theMaterial->getStrain();
}
