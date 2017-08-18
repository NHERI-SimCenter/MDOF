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
                                                                        
// $Revision: 1.14 $
// $Date: 2010-02-04 01:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLength.h,v $
                                                                        
                                                                        
#ifndef ZeroLength_h
#define ZeroLength_h

// Description: This file contains the class definition for ZeroLength.
// A ZeroLength element is defined by two nodes with the same coordinate.
// One or more material objects may be associated with the nodes to
// provide a force displacement relationship.
// ZeroLength element will work with 1d, 2d, or 3d material models.
//
// What: "@(#) ZeroLength.h, revA"

#include <Element.h>
#include <Matrix.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

// Type of dimension of element NxDy has dimension x=1,2,3 and
// y=2,4,6,12 degrees-of-freedom for the element
enum Etype { D1N2, D2N4, D2N6, D3N6, D3N12 };


class Node;
class Channel;
class UniaxialMaterial;
class Response;

class ZeroLength : public Element
{
  public:
    
  // Constructor for a single 1d material model
  ZeroLength(int tag, 			      
	     int Nd1, int Nd2, 
	     UniaxialMaterial& theMaterial,
	     double PdivL = 0.);

    ZeroLength();    
    ~ZeroLength();

    const char *getClassType(void) const {return "ZeroLength";};

    // public methods to obtain inforrmation about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);    

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    int setParameter(const char **argv, int argc, Parameter &param);
    

    // methods specific to MDOF
    double getForce();
    double getDrift();
  protected:
    
  private:
    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int numDOF;	                        // number of dof for ZeroLength
    Matrix transformation;		// transformation matrix for orientation
    Node *theNodes[2];

    Matrix *theMatrix; 	    	// pointer to objects matrix (a class Matrix)
    Vector *theVector;      	// pointer to objects vector (a class Vector)

    // Storage for uniaxial material models
    UniaxialMaterial *theMaterial;      // array of pointers to 1d materials

    // static data - single copy for all objects of the class	
    static Matrix ZeroLengthM2;   // class wide matrix for 2*2
    static Vector ZeroLengthV2;   // class wide Vector for size 2

    double PdivL;
};

#endif




