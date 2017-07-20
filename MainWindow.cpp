/* *****************************************************************************
Copyright (c) 2016-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: fmckenna

#include "MainWindow.h"
#include "ui_MainWindow.h"
#include "EarthquakeRecord.h"

#include <QDebug>
#include <QSlider>
//#include <QtNetwork>
#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QLabel>
#include <QPushButton>
#include <qcustomplot.h>
#include <MyGlWidget.h>


//#include <../widgets/InputSheetBM/SimpleSpreadsheetWidget.h>
#include <SimpleSpreadsheetWidget.h>
#include <NodeResponseWidget.h>

#include <QtNetwork/QNetworkAccessManager>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QNetworkRequest>

// OpenSees include files
#include <Node.h>
#include <ID.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <Domain.h>
#include <StandardStream.h>
//#include <LinearCrdTransf3d.h>
#include <Steel01.h>
#include <ZeroLength.h>
#include <LegendreBeamIntegration.h>
#include <ElasticSection3d.h>
#include <LinearSeries.h>
#include <NodalLoad.h>
#include <LoadPattern.h>
#include <SimulationInformation.h>
#include <PathSeries.h>
#include <GroundMotion.h>
#include <UniformExcitation.h>

#include <Newmark.h>
#include <RCM.h>
#include <PlainNumberer.h>
#include <NewtonRaphson.h>
#include <CTestNormDispIncr.h>
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSOE.h>
#ifdef _FORTRAN_LIBS
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#endif
#include <DirectIntegrationAnalysis.h>
#include <AnalysisModel.h>
#include "elCentro.AT2"
#include <Vector.h>
#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>

//style inludes
#include <QGroupBox>
#include <QFrame>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
Domain theDomain;

//
// procedure to create a QLabel QLineInput pair, returns pointer to QLineEdit created
//

QLineEdit *
createTextEntry(QString text,
                QVBoxLayout *theLayout,
                int minL=100,
                int maxL=100)
{
    QHBoxLayout *entryLayout = new QHBoxLayout();
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(text);

    QLineEdit *res = new QLineEdit();
    res->setMinimumWidth(minL);
    res->setMaximumWidth(maxL);
    res->setValidator(new QDoubleValidator);

    entryLayout->addWidget(entryLabel);
    entryLayout->addWidget(res);

    entryLayout->setSpacing(10);
    entryLayout->setMargin(0);

    theLayout->addLayout(entryLayout);


    return res;
}


//
// procedure to create a QLabel QLabel pair, returns pointer to QLabel created
//


QLabel *
createLabelEntry(QString text,
                 QVBoxLayout *theLayout,
                 int minL=100,
                 int maxL=100)
{
    QHBoxLayout *entryLayout = new QHBoxLayout();
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(text);

    QLabel *res = new QLabel();
    res->setMinimumWidth(minL);
    res->setMaximumWidth(maxL);

    entryLayout->addWidget(entryLabel);
    entryLayout->addWidget(res);

    entryLayout->setSpacing(10);
    entryLayout->setMargin(0);

    theLayout->addLayout(entryLayout);

    return res;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    numFloors(0), periods(0), buildingW(0), buildingH(1), storyK(0),
    weights(0), k(0), fy(0), b(0), dampRatios(0), floorHeights(0), storyHeights(0),
    dampingRatio(0.02), g(386.4), dt(0), gMotion(0),
    includePDelta(true), needAnalysis(true), eqData(0), dispResponses(0), maxDisp(1),
    movingSlider(false), fMinSelected(-1),fMaxSelected(-1), sMinSelected(-1),sMaxSelected(-1),
    time(1561),excitationValues(1561), graph(0), groupTracer(0),floorSelected(-1),storySelected(-1)
{

    createActions();

    // create a main layout
    mainLayout = new QHBoxLayout();
    largeLayout = new QVBoxLayout();

    // create input and output panel layouts and each to main layout
    createHeaderBox();
    createInputPanel();
    createOutputPanel();
    //createFooterBox();


    // create a widget in which to show everything //ALSO SET TO LARGE LAYOUT
    QWidget *widget = new QWidget();
    widget->setLayout(largeLayout);
    this->setCentralWidget(widget);


    //    resize(QDesktopWidget().availableGeometry(this).size() * 0.7);


    QRect rec = QApplication::desktop()->screenGeometry();

    int height = 0.7*rec.height();
    int width = 0.7*rec.width();

    this->resize(width, height);
    //
    // create 2 blank motions & make elCentro current
    //

    QStringList elCentrolist = elCentroTextData.split(QRegExp("[\r\n\t ]+"), QString::SkipEmptyParts);
    Vector *elCentroData = new Vector(elCentrolist.size()+1);
    // qDebug() << elCentrolist.size();

    (*elCentroData)(0) = 0;
    time[0]=0.;
    excitationValues[0]=0.;
    double maxValue = 0;
    for (int i = 0; i < elCentrolist.size(); ++i) {
        double value = elCentrolist.at(i).toDouble();
        (*elCentroData)(i+1) = value;
        time[i+1]=i*0.02;
        excitationValues[i+1]=value;
        if (fabs(value) > maxValue)
            maxValue = fabs(value);
    }
    dt = 0.02;
    numSteps = 1560;
    QString elCentroString("elCentro");
    EarthquakeRecord *elCentro = new EarthquakeRecord(elCentroString, 1560, 0.02, elCentroData);
    records.insert(std::make_pair(QString("elCentro"), elCentro));

    // create blank motion
    Vector *blankData = new Vector(100);
    QString blankString("Blank");
    EarthquakeRecord *blank = new EarthquakeRecord(blankString, 100, 0.02, blankData);
    records.insert(std::make_pair(blankString, blank));

    inMotion->addItem(elCentroString);
    inMotion->addItem(tr("Blank"));

    // create a basic model with defaults
    this->setBasicModel(5, 5*100, 5*144, 31.54, .05, 386.4);
    //setBasicModel(4,4,4,4,.02,386.4);

    // access a web page which will increment the usage count for this tool
    manager = new QNetworkAccessManager(this);

    connect(manager, SIGNAL(finished(QNetworkReply*)),
            this, SLOT(replyFinished(QNetworkReply*)));

    manager->get(QNetworkRequest(QUrl("http://opensees.berkeley.edu/OpenSees/developer/mdofUse.php")));
}

MainWindow::~MainWindow()
{
    delete ui;

    if (weights != 0)
        delete [] weights;
    if (k != 0)
        delete [] k;
    if (fy != 0)
        delete [] fy;
    if (b != 0)
        delete [] b;
    if (gMotion != 0)
        delete [] gMotion;
    if (storyHeights != 0)
        delete [] storyHeights;
    if (floorHeights != 0)
        delete [] floorHeights;
    if (dampRatios != 0)
        delete [] dampRatios;
}

void MainWindow::draw(MyGlWidget *theGL)
{
    if (needAnalysis == true) {
        doAnalysis();
    }
    theGL->reset();

    for (int i=0; i<numFloors; i++) {
        if (i == storySelected)
            theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 0, 0, 0);
        else if (i >= sMinSelected && i <= sMaxSelected)
            theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 1, 0, 0);
        else
            theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 0, 0, 0);
    }
    
    for (int i=0; i<=numFloors; i++) {
        if (i == floorSelected)
            theGL->drawPoint(i, dispResponses[i][currentStep],floorHeights[i], 10, 0, 0, 0);
        else if (i >= fMinSelected && i <= fMaxSelected)
            theGL->drawPoint(i, dispResponses[i][currentStep],floorHeights[i], 10, 1, 0, 0);
        else
            theGL->drawPoint(i, dispResponses[i][currentStep],floorHeights[i], 10, 0, 0, 1);
    }

    // display range of displacement
    static char maxDispString[30];
    snprintf(maxDispString, 50, "%.3e", maxDisp);
    theGL->drawLine(0, -maxDisp, 0.0, maxDisp, 0.0, 1.0, 0., 0., 0.);

    currentTime->setText(QString().setNum(currentStep*dt,'f',2));
    currentDisp->setText(QString().setNum(dispResponses[numFloors][currentStep],'f',2));
    theGL->drawBuffers();

    // update red dot on earthquake plot
    /*
    groupTracer->setGraph(0);
    groupTracer->setGraph(graph);
    groupTracer->setGraphKey(currentStep*dt);
    groupTracer->updatePosition();
    //earthquakePlot->replot();

*/
}



void MainWindow::updatePeriod()
{
    //periods = 2.0;
}

void MainWindow::setBasicModel(int numF, double W, double H, double K, double zeta, double grav)
{
    if (numFloors != numF) {

        // if invalid numFloor, return
        if (numF <= 0)
            return;

        // resize arrays
        if (weights != 0)
            delete [] weights;
        if (k != 0)
            delete [] k;
        if (fy != 0)
            delete [] fy;
        if (b != 0)
            delete [] b;
        if (floorHeights != 0)
            delete [] floorHeights;
        if (storyHeights != 0)
            delete [] storyHeights;
        if (dampRatios != 0)
            delete [] dampRatios;

        if (dispResponses != 0) {
            for (int j=0; j<numFloors+1; j++)
                delete [] dispResponses[j];
            delete [] dispResponses;
        }

        weights = new double[numF];
        k = new double[numF];
        fy = new double[numF];
        b = new double[numF];
        floorHeights = new double[numF+1];
        storyHeights = new double[numF];
        dampRatios = new double[numF];

        dispResponses = new double *[numF+1];

        //for (int i=0; i<numF+1; i++) {}
        //numSteps = 2000;
        for (int i=0; i<numF+1; i++) {
            dispResponses[i] = new double[numSteps+1]; // +1 as doing 0 at start
        }
    }

    // set values
    double floorW = W/(numF);
    buildingW = W;
    storyK = K;
    numFloors = numF;
    buildingH = H;
    dampingRatio = zeta;
    g=grav;

    for (int i=0; i<numF; i++) {
        weights[i] = floorW;
        k[i] = K;
        fy[i] = 1.0e100;
        b[i] = 0.;
        floorHeights[i] = i*buildingH/(1.*numF);
        storyHeights[i] = buildingH/(1.*numF);
        dampRatios[i] = zeta;
    }
    floorHeights[numF] = buildingH;

    this->updatePeriod();

    // update text boxes
    //inPeriod->setText(QString::number(period));

    inWeight->setText(QString::number(buildingW));
    inK->setText(QString::number(storyK));
    inFloors->setText(QString::number(numF));
    inHeight->setText(QString::number(buildingH));
    inDamping->setText(QString::number(zeta));
    inGravity->setText(QString::number(g));
    needAnalysis = true;
    this->reset();

    theNodeResponse->setFloor(numF);
}


void
MainWindow::on_includePDeltaChanged(int state)
{
    if (state == Qt::Checked)
        includePDelta = true;
    else
        includePDelta = false;

    this->reset();
}

void MainWindow::on_inFloors_editingFinished()
{
    QString textFloors =  inFloors->text();
    int numFloorsText = textFloors.toInt();
    if (numFloorsText != numFloors) {
        this->setBasicModel(numFloorsText, buildingW, numFloorsText, storyK, dampingRatio, 386.4);
    }
    //  inWeight->setFocus();

    this->reset();
}

void MainWindow::on_inWeight_editingFinished()
{
    QString textW =  inWeight->text();
    double textToDoubleW = textW.toDouble();
    if (textToDoubleW != buildingW) {
        // set values
        double floorW = textToDoubleW/(numFloors);
        for (int i=0; i<numFloors; i++) {
            weights[i] = floorW;
        }
    }
    // inHeight->setFocus();
    this->reset();
}

void MainWindow::on_inHeight_editingFinished()
{
    QString textH =  inHeight->text();
    double textToDoubleH = textH.toDouble();
    if (textToDoubleH != buildingH) {
        // set values
        buildingH = textToDoubleH;
        double deltaH = buildingH/numFloors;
        floorHeights[0] = 0;
        for (int i=0; i<numFloors; i++) {
            storyHeights[i] = deltaH;
            floorHeights[i+1] = deltaH + floorHeights[i];
        }
    }
    //   inK->setFocus();
    this->reset();
}


void MainWindow::on_inK_editingFinished()
{
    QString text =  inK->text();
    double textToDouble = text.toDouble();
    for (int i=0; i<numFloors; i++)
        k[i] = textToDouble;
    //  inDamping->setFocus();
    this->reset();
}

void MainWindow::on_inDamping_editingFinished()
{
    QString text =  inDamping->text();
    double textToDouble = text.toDouble();
    dampingRatio = textToDouble;
    //    inGravity->setFocus();
    for (int i=0; i<numFloors; i++) {
        dampRatios[i]=dampingRatio;
    }
    this->reset();
}

void MainWindow::on_inFloorWeight_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inFloorWeight->text();
    double textToDouble = text.toDouble();
    for (int i=fMinSelected; i<=fMaxSelected; i++)
        weights[i] = textToDouble;

    buildingW = 0;
    for (int i=0; i<numFloors; i++)
        buildingW = buildingW+weights[i];

    inWeight->setText(QString::number(buildingW));
    this->reset();
}




void MainWindow::on_inGravity_editingFinished()
{
    QString text =  inGravity->text();
    double textToDouble = text.toDouble();
    g = textToDouble;

    this->reset();
}

void MainWindow::on_inStoryHeight_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryHeight->text();
    double newStoryHeight = text.toDouble();
    double currentStoryHeight = 0;
    double *newFloorHeights = new double[numFloors+1];

    // determine new floor heights, cludgy can rewrite now store storyHeights
    newFloorHeights[0] = 0;
    for (int i=0; i<sMinSelected; i++)
        newFloorHeights[i+1] = floorHeights[i+1];

    for (int i=sMinSelected; i<=sMaxSelected; i++) {
        newFloorHeights[i+1] = newFloorHeights[i]+newStoryHeight;
        storyHeights[i] = newStoryHeight;
    }

    for (int i=sMaxSelected+1; i<numFloors; i++)
        newFloorHeights[i+1] = newFloorHeights[i]+floorHeights[i+1]-floorHeights[i];

    // delete old array and reset pointer
    buildingH = newFloorHeights[numFloors];
    inHeight->setText(QString::number(buildingH));

    delete [] floorHeights;
    floorHeights = newFloorHeights;

    // move focus, update graphic and set analysis flag
    inStoryK->setFocus();
    this->reset();
}

void MainWindow::on_inStoryK_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryK->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        k[i] = textToDouble;

    inStoryFy->setFocus();
    this->reset();
}

void MainWindow::on_inStoryFy_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryFy->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        fy[i] = textToDouble;

    inStoryB->setFocus();
    this->reset();
}

void MainWindow::on_inStoryB_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryB->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        b[i] = textToDouble;

    inStoryHeight->setFocus();
    this->reset();
}


void MainWindow::doAnalysis()
{

    if (needAnalysis == true) {

        // clear existinqDebugg model
        theDomain.clearAll();
        OPS_clearAllUniaxialMaterial();
        ops_Dt = 0.0;

        Node **theNodes = new Node *[numFloors+1];
        for (int i=0; i<= numFloors; i++) {
            Matrix theMass(1,1);
            Node *theNode=new Node(i+1,1, 0.0);
            theNodes[i] = theNode;
            theDomain.addNode(theNode);
            if (i == 0) {
                SP_Constraint *theSP = new SP_Constraint(i+1, 0, 0., true);
                theDomain.addSP_Constraint(theSP);
            } else {
                theMass(0,0) = weights[i-1]/g;
                theNode->setMass(theMass);
            }
        }


        // create the vectors for the element orientation
        Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
        Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

        double axialLoad = 0;
        for (int i=numFloors;  i>0; i--) {
            UniaxialMaterial *theMat = new Steel01(i,fy[i-1],k[i-1],b[i-1]);
            //  ZeroLength *theEle = new ZeroLength(i+1, 1, i+1, i+2,
            //x, y, *theMat, 0);
            double PdivL = 0.0;
            if (includePDelta == true && storyHeights[i-1] != 0) {
                axialLoad = axialLoad + weights[i-1];
                PdivL = -axialLoad/storyHeights[i-1]; // negative for compression
            }
            ZeroLength *theEle = new ZeroLength(i, i, i+1, *theMat, PdivL);
            theDomain.addElement(theEle);
        }

        //
        // create load pattern and add loads
        //

        PathSeries *theSeries = new PathSeries(1, *eqData, dt, g);
        GroundMotion *theGroundMotion = new GroundMotion(0,0,theSeries);
        LoadPattern *theLoadPattern = new UniformExcitation(*theGroundMotion, 0, 1);
        //   theLoadPattern->setTimeSeries(theTimeSeries);
        //   static Vector load(1); load.Zero(); load(0) = 1;
        //   NodalLoad *theLoad = new NodalLoad(0, numFloors, load);
        //   theLoadPattern->addNodalLoad(theLoad);
        theDomain.addLoadPattern(theLoadPattern);

        //theDomain.Print(opserr);
        //
        // create the analysis
        //

        AnalysisModel     *theModel = new AnalysisModel();
        CTestNormDispIncr *theTest = new CTestNormDispIncr(1.0e-3, 20, 0);
        EquiSolnAlgo      *theSolnAlgo = new NewtonRaphson(INITIAL_TANGENT);
        TransientIntegrator  *theIntegrator = new Newmark(0.5, 0.25);
        //ConstraintHandler *theHandler = new TransformationConstraintHandler();
        ConstraintHandler *theHandler = new PlainHandler();
        RCM               *theRCM = new RCM();
        DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);
#ifdef _FORTRAN_LIBS
        BandGenLinSolver  *theSolver = new BandGenLinLapackSolver();
        LinearSOE         *theSOE = new BandGenLinSOE(*theSolver);
#else
        ProfileSPDLinSolver  *theSolver = new ProfileSPDLinDirectSolver();
        LinearSOE         *theSOE = new ProfileSPDLinSOE(*theSolver);
#endif
        DirectIntegrationAnalysis    theAnalysis(theDomain,
                                                 *theHandler,
                                                 *theNumberer,
                                                 *theModel,
                                                 *theSolnAlgo,
                                                 *theSOE,
                                                 *theIntegrator);
        theSolnAlgo->setConvergenceTest(theTest);

        SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
        EigenSOE *theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theModel);
        theAnalysis.setEigenSOE(*theEigenSOE);

        int ok = theAnalysis.eigen(numFloors,true);
        const Vector &theEig = theDomain.getEigenvalues();
        if (ok == 0)
          for (int i=0; i<numFloors; i++)
            if (theEig(i) <= 0)
                ok = -1;

        if (ok != 0) {
            QMessageBox::warning(this, tr("Application"),
                                 tr("Eigenvalue Analysis Failed. Possible Causes: Negative stiffness "
                                    "(either due to negative story stiffness value or large Axial force leading to "
                                    "large negative PDelta contribibution"));
            needAnalysis = false;
                                // .arg(QDir::toNativeSeparators(fileName), file.errorString()));
            return;
        }


        Vector dampValues(numFloors);
        for (int i=0; i<numFloors; i++) {
            dampValues(i)=dampRatios[i];
        }
        theDomain.setModalDampingFactors(&dampValues);


        double T1 = 2*3.14159/sqrt(theEig(0));
        //qDebug() << T1;
        //inPeriod->setText(QString::number(T1));



        //
        //analyze & get results
        //
        maxDisp = 0;
        for (int i=0; i<=numSteps; i++) { // <= due to adding 0 at start
            int ok = theAnalysis.analyze(1, dt);
            if (ok != 0) {
                QMessageBox::warning(this, tr("Application"),
                                     tr("Transient Analysis Failed"));
                needAnalysis = false;

                                    // .arg(QDir::toNativeSeparators(fileName), file.errorString()));
                break;
            }
            for (int j=0; j<numFloors+1; j++) {
                double nodeDisp = theNodes[j]->getDisp()(0);
                dispResponses[j][i] = nodeDisp;
                if (fabs(nodeDisp) > maxDisp)
                    maxDisp = fabs(nodeDisp);
            }
            if (ok != 0)
                break;
        }

        // clean up memory
        delete [] theNodes;
        maxDispLabel->setText(QString().setNum(maxDisp,'f',2));
        currentPeriod->setText(QString().setNum(T1,'f',2));
        // reset values, i.e. slider position, current displayed step, and display properties
        needAnalysis = false;
        currentStep = 0;
        //  groupTracer->setGraphKey(0);
        slider->setSliderPosition(0);
        myGL->update();

        int nodeResponseFloor = theNodeResponse->getFloor();
        nodeResponseValues.resize(numSteps);

        for (int i = 0; i < numSteps; ++i) {
            nodeResponseValues[i]=dispResponses[nodeResponseFloor][i];
        }
        theNodeResponse->setData(nodeResponseValues,time,numSteps,dt);
    }
}

void
MainWindow::setFloorResponse(int floor)
{
    if (floor > 0 && floor <= numFloors) {
        for (int i = 0; i < numSteps; ++i) {
            nodeResponseValues[i]=dispResponses[floor][i];
        }
        theNodeResponse->setData(nodeResponseValues,time,numSteps,dt);
    }
}

void MainWindow::reset() {

    needAnalysis = true;
    myGL->update();

    // update the properties table

    theSpreadsheet->clear();
    theSpreadsheet->setColumnCount(6);
    theSpreadsheet->setRowCount(numFloors);
    theSpreadsheet->horizontalHeader()->setStretchLastSection(true);// horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
    //theSpreadsheet->setFixedWidth(344);
    updatingPropertiesTable = true;
    theSpreadsheet->setHorizontalHeaderLabels(QString(" Weight ; Height ;    K    ;    Fy    ;    b    ;  zeta").split(";"));
    for (int i=0; i<numFloors; i++) {
        theSpreadsheet->setItem(i,0,new QTableWidgetItem(QString().setNum(weights[i])));
        theSpreadsheet->setItem(i,1,new QTableWidgetItem(QString().setNum(storyHeights[i])));
        theSpreadsheet->setItem(i,2,new QTableWidgetItem(QString().setNum(k[i])));
        theSpreadsheet->setItem(i,3,new QTableWidgetItem(QString().setNum(fy[i])));
        theSpreadsheet->setItem(i,4,new QTableWidgetItem(QString().setNum(b[i])));
        theSpreadsheet->setItem(i,5,new QTableWidgetItem(QString().setNum(dampRatios[i])));
    }
    theSpreadsheet->resizeRowsToContents();
    theSpreadsheet->resizeColumnsToContents();

    updatingPropertiesTable = false;

    floorSelected = -1;
    storySelected = -1;
}

void MainWindow::on_theSpreadsheet_cellChanged(int row, int column)
{
    if (updatingPropertiesTable == false) {
        QString text = theSpreadsheet->item(row,column)->text();

        bool ok;
        double textToDouble = text.toDouble(&ok);
        if (column == 0) {
            weights[row] = textToDouble;
            buildingW = 0;
            for (int i=0; i<numFloors; i++)
                buildingW += weights[i];
            inWeight->setText(QString::number(buildingW));
        } else if (column  == 1) {
            storyHeights[row] = textToDouble;
            for (int i=row; i<numFloors; i++)
                floorHeights[i+1] = floorHeights[i]+storyHeights[i];
            buildingH = floorHeights[numFloors];
            inHeight->setText(QString::number(buildingH));
        } else if (column == 2) {
            k[row] = textToDouble;
        } else if (column == 3) {
            fy[row] = textToDouble;
        } else if (column == 4) {
            b[row] = textToDouble;
        } else
            dampRatios[row] = textToDouble;


        needAnalysis = true;
        myGL->update();
    }
}


void MainWindow::on_stopButton_clicked()
{
    stopRun = true;
}

void MainWindow::on_exitButton_clicked()
{
    close();
}

void MainWindow::on_runButton_clicked()
{
    stopRun = false;
    if (needAnalysis == true) {
        this->doAnalysis();
    }

    //currentStep = 0;
    do { //while (currentStep < numSteps && stopRun == false){

        slider->setSliderPosition(currentStep);
        myGL->repaint();
        QCoreApplication::processEvents();

        currentStep++;

    } while (currentStep <= numSteps && stopRun == false); // <= added 0 to ground motion
}


void MainWindow::on_slider_valueChanged(int value)
{
    if (movingSlider == true) {
        stopRun = true;
        if (needAnalysis == true) {
            this->doAnalysis();
            myGL->update();
        }
        currentStep = slider->value();

        myGL->repaint();
    }
}

void MainWindow::on_slider_sliderPressed()
{
    movingSlider = true;
}

void MainWindow::on_slider_sliderReleased()
{
    movingSlider = false;
}

float
MainWindow::setSelectionBoundary(float y1, float y2)
{
    // determine min and max of y1,y2
    float yMin = 0;
    float yMax = 0;
    if (y1 < y2) {
        yMin = y1;
        yMax = y2;
    } else {
        yMin = y2;
        yMax = y1;
    }

    // determine min and max of nodes in [y1,y2] range
    fMinSelected = -1;
    fMaxSelected = -1;

    for (int i=0; i<=numFloors; i++) {
        if (floorHeights[i] >= yMin && floorHeights[i] <= yMax) {
            if (fMinSelected == -1)
                fMinSelected = i;
            fMaxSelected = i;
        }
    }

    // determine min and max of stories in [y1, y2] range;
    sMinSelected = -1;
    sMaxSelected = -1;
    for (int i=0; i<numFloors; i++) {
        double midStoryHeight = (floorHeights[i]+floorHeights[i+1])/2.;
        if (midStoryHeight >= yMin && midStoryHeight <= yMax) {
            if (sMinSelected == -1)
                sMinSelected = i;
            sMaxSelected = i;
        }
    }
    //   qDebug() << "sMinSelected: " << sMinSelected << " sMaxSelected: " << sMaxSelected;
    //   qDebug() << "fMinSelected: " << fMinSelected << " fMaxSelected: " << fMaxSelected;

    updatingPropertiesTable = true;

    if (fMinSelected == 0 && fMaxSelected == numFloors) {

        floorMassFrame->setVisible(false);
        storyPropertiesFrame->setVisible(false);
        spreadSheetFrame->setVisible(true);

    } else if (fMinSelected == fMaxSelected && fMinSelected != -1) {

        floorMassFrame->setVisible(true);
        if (sMinSelected == -1 && sMaxSelected == -1)
            storyPropertiesFrame->setVisible(false);
        else
            storyPropertiesFrame->setVisible(true);
        spreadSheetFrame->setVisible(false);
        floorSelected=-1;
        storySelected=-1;

    } else if (fMinSelected != -1 && fMaxSelected != -1) {
        floorMassFrame->setVisible(true);
        storyPropertiesFrame->setVisible(true);
        spreadSheetFrame->setVisible(false);
        floorSelected=-1;
        storySelected=-1;

    } else if (sMinSelected != -1 && sMaxSelected != -1) {
        floorMassFrame->setVisible(false);
        storyPropertiesFrame->setVisible(true);
        spreadSheetFrame->setVisible(false);
        floorSelected=-1;
        storySelected=-1;

    } else {

        floorMassFrame->setVisible(false);
        storyPropertiesFrame->setVisible(false);
        spreadSheetFrame->setVisible(false);
        floorSelected=-1;
        storySelected=-1;
    }

    updatingPropertiesTable = false;

    //
    // based on min, max nodes enable/disable lineEdits & set text
    //

    myGL->repaint();
    return 0;
}




void MainWindow::on_theSpreadsheet_cellClicked(int row, int column)
{
    if (column == 0) {
        floorSelected = row+1;
        storySelected = -1;
    }
    else if (column > 0 && column < 5) {
        storySelected = row;
        floorSelected = -1;
    } else {
        storySelected = -1;
        floorSelected = -1;
    }

    myGL->repaint();
}

void MainWindow::replyFinished(QNetworkReply *pReply)
{
    qDebug() << "INTERNET HIT";

    //QByteArray data=pReply->readAll();
    //QString str(data);
}



void MainWindow::on_inMotionSelection_currentTextChanged(const QString &arg1)
{
    std::map<QString, EarthquakeRecord *>::iterator it;
    it = records.find(arg1);
    if (it != records.end()) {
        EarthquakeRecord *theRecord = records.at(arg1);
        numSteps =  theRecord->numSteps;
        dt =  theRecord->dt;
        eqData = theRecord->data;
        double maxValue = 0;
        excitationValues.resize(numSteps);
        time.resize(numSteps);
        for (int i = 0; i < numSteps; ++i) {
            double value = (*eqData)[i];
            time[i]=i*dt;
            excitationValues[i]=value;
            if (fabs(value) > maxValue)
                maxValue = fabs(value);
        }

        // reset earthquake plot
        earthquakePlot->clearGraphs();
        graph = earthquakePlot->addGraph();
        earthquakePlot->graph(0)->setData(time, excitationValues);
        earthquakePlot->xAxis->setRange(0, numSteps*dt);
        earthquakePlot->yAxis->setRange(-maxValue, maxValue);
        earthquakePlot->axisRect()->setAutoMargins(QCP::msNone);
        earthquakePlot->axisRect()->setMargins(QMargins(0,0,0,0));
        /*
        if (groupTracer != 0)
            delete groupTracer;
        groupTracer = new QCPItemTracer(earthquakePlot);
        groupTracer->setGraph(graph);
        groupTracer->setGraphKey(0);
        groupTracer->setInterpolating(true);
        groupTracer->setStyle(QCPItemTracer::tsCircle);
        groupTracer->setPen(QPen(Qt::red));
        groupTracer->setBrush(Qt::red);
        groupTracer->setSize(7);
*/
        // reset slider range
        slider->setRange(0, numSteps);

        this->reset();
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    QApplication::quit();
}


bool MainWindow::save()
{
    if (currentFile.isEmpty()) {
        return saveAs();
    } else {
        return saveFile(currentFile);
    }
}

bool MainWindow::saveAs()
{
    //
    // get filename
    //

    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if (dialog.exec() != QDialog::Accepted)
        return false;

    // and save the file
    return saveFile(dialog.selectedFiles().first());
}

void MainWindow::open()
{
    QString fileName = QFileDialog::getOpenFileName(this);
    if (!fileName.isEmpty())
        loadFile(fileName);
}

void MainWindow::newFile()
{
    // clear old
    //inputWidget->clear();

    // set currentFile blank
    setCurrentFile(QString());
}


void MainWindow::setCurrentFile(const QString &fileName)
{
    currentFile = fileName;
    //  setWindowModified(false);

    QString shownName = currentFile;
    if (currentFile.isEmpty())
        shownName = "untitled.json";

    setWindowFilePath(shownName);
}


void MainWindow::on_addMotion_clicked()
{

    //
    // open files
    //

    QString inputMotionName = QFileDialog::getOpenFileName(this);
    if (inputMotionName.isEmpty())
        return;

    QFile file(inputMotionName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(inputMotionName), file.errorString()));
        return;
    }

    // place contents of file into json object
    QString val;
    val=file.readAll();
    QJsonDocument doc = QJsonDocument::fromJson(val.toUtf8());

    QJsonObject jsonObject = doc.object();

    QJsonValue theValue = jsonObject["name"];
    QString name =theValue.toString();
    theValue = jsonObject["numPoints"];
    int numPoints = theValue.toInt();
    theValue = jsonObject["dT"];
    double dT = theValue.toDouble();
    theValue = jsonObject["data"];
    QJsonArray data = theValue.toArray();

    // create blank motion
    Vector *theData = new Vector(numPoints);

    for (int i=0; i<numPoints; i++) {
        theValue = data.at(i);
        (*theData)[i] = theValue.toDouble();
    }

    EarthquakeRecord *theRecord = new EarthquakeRecord(name, numPoints, dT, theData);
    records.insert(std::make_pair(name, theRecord));

    inMotion->addItem(name);

    // close file
    file.close();
}

bool MainWindow::saveFile(const QString &fileName)
{
    //
    // open file
    //

    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName),
                                  file.errorString()));
        return false;
    }


    //
    // create a json object, fill it in & then use a QJsonDocument
    // to write the contents of the object to the file in JSON format
    //

    QJsonObject json;
    json["numFloors"]=numFloors;
    json["buildingHeight"]=buildingH;
    json["buildingWeight"]=buildingW;
    json["K"]=storyK;
    json["dampingRatio"]=dampingRatio;
    json["G"]=g;
    json["currentMotion"]=inMotion->currentText();
    json["currentMotionIndex"]=inMotion->currentIndex();

    QJsonArray weightsArray;
    QJsonArray kArray;
    QJsonArray fyArray;
    QJsonArray bArray;
    QJsonArray heightsArray;
    QJsonArray dampArray;

    for (int i=0; i<numFloors; i++) {
        weightsArray.append(weights[i]);
        kArray.append(k[i]);
        fyArray.append(fy[i]);
        bArray.append(b[i]);
        heightsArray.append(storyHeights[i]);
        dampArray.append(dampRatios[i]);
    }


    json["floorWeights"]=weightsArray;
    json["storyK"]=kArray;
    json["storyFy"]=fyArray;
    json["storyB"]=bArray;
    json["storyHeights"]=heightsArray;
    json["dampRatios"]=dampArray;

    //inputWidget->outputToJSON(json);
    QJsonDocument doc(json);
    file.write(doc.toJson());

    // close file
    file.close();

    // set current file
    setCurrentFile(fileName);

    return true;
}

void MainWindow::copyright()
{
    QString textCopyright = "\
            Copyright (c) 2017-2018, The Regents of the University of California (Regents).\
            All rights reserved.\
            <p>\
            Redistribution and use in source and binary forms, with or without \
            modification, are permitted provided that the following conditions are met:\
            <p>\
            1. Redistributions of source code must retain the above copyright notice, this\
               list of conditions and the following disclaimer.\
            2. Redistributions in binary form must reproduce the above copyright notice,\
               this list of conditions and the following disclaimer in the documentation\
               and/or other materials provided with the distribution.\
            <p>\
            THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND\
            ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\
            WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\
            DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\
            ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\
            (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\
            LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\
            ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\
            (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\
            SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\
            <p>\
            The views and conclusions contained in the software and documentation are those\
            of the authors and should not be interpreted as representing official policies,\
            either expressed or implied, of the FreeBSD Project.\
            <p>\
            REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, \
            THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.\
            THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS \
            PROVIDED \"AS IS\". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,\
            UPDATES, ENHANCEMENTS, OR MODIFICATIONS.\
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            This makes use of the QT packages (unmodified): core, gui, widgets and network\
            <p>\
            QT is copyright \"The Qt Company Ltd&quot; and licensed under the GNU Lesser General \
            Public License (version 3) which references the GNU General Public License (version 3)\
            <p>\
            These Licenses can be found at: &lt;http://www.gnu.org/licenses/&gt;";


   QMessageBox msgBox;
   QSpacerItem *theSpacer = new QSpacerItem(500, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
   msgBox.setText(textCopyright);
   QGridLayout *layout = (QGridLayout*)msgBox.layout();
   layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
   msgBox.exec();

}


void MainWindow::version()
{
   QMessageBox::about(this, tr("Version"),
            tr("Version 0.1 Beta Release "));
}

void MainWindow::about()
{
    QString textAbout = "\
   This is the Multiple Degree of Freedom (MDOF) tool\
   It presents a shear spring model of a multi-story building\
   All units are in sec, kips, inches.\
   <p>\
   This tool is in beta release mode\
   Any suggestions or bugs should be submitted to \
   &lthttps://github.com/NHERI-SimCenter/MDOF/issues&gt\
   ";

   QMessageBox msgBox;
   QSpacerItem *theSpacer = new QSpacerItem(500, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
   msgBox.setText(textAbout);
   QGridLayout *layout = (QGridLayout*)msgBox.layout();
   layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
   msgBox.exec();
}

void MainWindow::submitFeedback()
{
   QDesktopServices::openUrl(QUrl("https://github.com/NHERI-SimCenter/MDOF/issues", QUrl::TolerantMode));
}

void MainWindow::loadFile(const QString &fileName)
{
    //
    // open files
    //

    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName), file.errorString()));
        return;
    }

    // place contents of file into json object
    QString val;
    val=file.readAll();
    QJsonDocument doc = QJsonDocument::fromJson(val.toUtf8());
    QJsonObject jsonObject = doc.object();

    QJsonValue theValue = jsonObject["numFloors"];
    numFloors=theValue.toInt();
    inFloors->setText(QString::number(numFloors));

    theValue = jsonObject["buildingHeight"];
    buildingH=theValue.toDouble();
    inHeight->setText(QString::number(buildingH));

    theValue = jsonObject["buildingWeight"];
    buildingW=theValue.toDouble();
    inWeight->setText(QString::number(buildingW));

    theValue = jsonObject["K"];
    storyK=theValue.toDouble();
    inK->setText(QString::number(storyK));

    theValue = jsonObject["G"];
    g=theValue.toDouble();
    inGravity->setText(QString::number(g));

    theValue = jsonObject["dampingRatio"];
    dampingRatio=theValue.toDouble();
    inDamping->setText(QString::number(dampingRatio));

    theValue = jsonObject["currentMotionIndex"];
    inMotion->setCurrentIndex(theValue.toInt());
    theValue = jsonObject["currentMotion"];

    if (weights != 0)
        delete [] weights;
    if (k != 0)
        delete [] k;
    if (fy != 0)
        delete [] fy;
    if (b != 0)
        delete [] b;
    if (floorHeights != 0)
        delete [] floorHeights;
    if (storyHeights != 0)
        delete [] storyHeights;
    if (dampRatios != 0)
        delete [] dampRatios;

    weights = new double[numFloors];
    k = new double[numFloors];
    fy = new double[numFloors];
    b = new double[numFloors];
    floorHeights = new double[numFloors+1];
    storyHeights = new double[numFloors];
    dampRatios = new double[numFloors];

    QJsonArray theArray;

    theValue = jsonObject["floorWeights"];
    theArray=theValue.toArray();

    for (int i=0; i<numFloors; i++)
        weights[i] = theArray.at(i).toDouble();

    theValue = jsonObject["storyK"];
    theArray=theValue.toArray();

    for (int i=0; i<numFloors; i++)
        k[i] = theArray.at(i).toDouble();

    theValue = jsonObject["storyFy"];
    theArray=theValue.toArray();

    for (int i=0; i<numFloors; i++)
        fy[i] = theArray.at(i).toDouble();

    theValue = jsonObject["storyB"];
    theArray=theValue.toArray();

    for (int i=0; i<numFloors; i++)
        b[i] = theArray.at(i).toDouble();

    theValue = jsonObject["storyHeights"];
    theArray=theValue.toArray();

    floorHeights[0] = 0;
    for (int i=0; i<numFloors; i++) {
        storyHeights[i] = theArray.at(i).toDouble();
        floorHeights[i+1] = floorHeights[i] + storyHeights[i];
    }

    theValue = jsonObject["dampRatios"];
    theArray=theValue.toArray();

    for (int i=0; i<numFloors; i++)
        dampRatios[i] = theArray.at(i).toDouble();


    dispResponses = new double *[numFloors+1];
    for (int i=0; i<numFloors+1; i++) {
        dispResponses[i] = new double[numSteps+1]; // +1 as doing 0 at start
    }


    this->reset();
    this->on_inMotionSelection_currentTextChanged(theValue.toString());
    theNodeResponse->setFloor(numFloors);

    // close file
    file.close();

    // given the json object, create the C++ objects
    // inputWidget->inputFromJSON(jsonObj);

    //setCurrentFile(fileName);
}


void MainWindow::createActions() {

    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));

    //const QIcon openIcon = QIcon::fromTheme("document-open", QIcon(":/images/open.png"));
    //const QIcon saveIcon = QIcon::fromTheme("document-save", QIcon(":/images/save.png"));

    //QToolBar *fileToolBar = addToolBar(tr("File"));

    QAction *newAction = new QAction(tr("&New"), this);
    newAction->setShortcuts(QKeySequence::New);
    newAction->setStatusTip(tr("Create a new file"));
    connect(newAction, &QAction::triggered, this, &MainWindow::newFile);
    fileMenu->addAction(newAction);
    //fileToolBar->addAction(newAction);


    QAction *openAction = new QAction(tr("&Open"), this);
    openAction->setShortcuts(QKeySequence::Open);
    openAction->setStatusTip(tr("Open an existing file"));
    connect(openAction, &QAction::triggered, this, &MainWindow::open);
    fileMenu->addAction(openAction);
    //fileToolBar->addAction(openAction);

    QAction *saveAction = new QAction(tr("&Save"), this);
    saveAction->setShortcuts(QKeySequence::Save);
    saveAction->setStatusTip(tr("Save the document to disk"));
    connect(saveAction, &QAction::triggered, this, &MainWindow::save);
    fileMenu->addAction(saveAction);

    QAction *saveAsAction = new QAction(tr("&Save As"), this);
    saveAction->setStatusTip(tr("Save the document with new filename to disk"));
    connect(saveAction, &QAction::triggered, this, &MainWindow::save);
    fileMenu->addAction(saveAsAction);

    // strangely, this does not appear in menu (at least on a mac)!! ..
    // does Qt not allow as in tool menu by default?
    // check for yourself by changing Quit to drivel and it works
    QAction *exitAction = new QAction(tr("&Quit"), this);
    connect(exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));
    // exitAction->setShortcuts(QKeySequence::Quit);
    exitAction->setStatusTip(tr("Exit the application"));
    fileMenu->addAction(exitAction);

    QMenu *viewMenu = menuBar()->addMenu(tr("&View"));

    theNodeResponse = new NodeResponseWidget(this);
    QDockWidget *nodeResponseDock = new QDockWidget(tr("Node Response"), this);
    //dockO->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea |\
    Qt::TopDockWidgetArea);
    nodeResponseDock->setWidget(theNodeResponse);
    nodeResponseDock->setAllowedAreas(Qt::NoDockWidgetArea);
    nodeResponseDock->setFloating(true);
    nodeResponseDock->close();

    viewMenu->addAction(nodeResponseDock->toggleViewAction());

   QMenu *helpMenu = menuBar()->addMenu(tr("&About"));
   QAction *infoAct = helpMenu->addAction(tr("&Information"), this, &MainWindow::about);
   QAction *submitAct = helpMenu->addAction(tr("&Provide Feedback"), this, &MainWindow::submitFeedback);
    //aboutAct->setStatusTip(tr("Show the application's About box"));
   QAction *aboutAct = helpMenu->addAction(tr("&Version"), this, &MainWindow::version);
    //aboutAct->setStatusTip(tr("Show the application's About box"));


   QAction *copyrightAct = helpMenu->addAction(tr("&Copyright"), this, &MainWindow::copyright);
    //aboutAct->setStatusTip(tr("Show the application's About box"));

}

void MainWindow::viewNodeResponse(){

}

void MainWindow::viewStoryResponse(){

}

void MainWindow::createHeaderBox() {

    //
    // Make the header layout
    // styleSheet

    headerLayout = new QHBoxLayout;

    QGroupBox *header =new QGroupBox(tr("Multiple Degrees of Freedom Application"));
    headerLayout->addWidget(header);

    largeLayout->addLayout(headerLayout);
}

void MainWindow::createFooterBox() {

    //
    // Make the footer layout
    // styleSheet

    QGroupBox *footer =new QGroupBox();
    QLabel *nsfLogo = new QLabel();
    QPixmap pixmap(":/mdof.gif");
    QPixmap newPixmap = pixmap.scaled(QSize(40,40),  Qt::KeepAspectRatio);
    nsfLogo->setPixmap(newPixmap);
    nsfLogo->setMask(newPixmap.mask());
    nsfLogo->show();

//    QLabel *simLogo = new QLabel();
//    QPixmap pixmap1("/Users/TylerDurden/Projects/sim/mdof_fork/simcenter_cut.png");
//    QPixmap simPixmap = pixmap1.scaled(QSize(40,40),  Qt::KeepAspectRatio);
//    simLogo->setPixmap(simPixmap);
//    simLogo->setMask(simPixmap.mask());
//    simLogo->show();

    QLabel *nsfText = new QLabel();
    nsfText->setObjectName(QString::fromUtf8("nsfText"));
    nsfText->setText(tr("This work is based on material supported by the National Science Foundation under grant 1612843-2"));

    footerLayout = new QHBoxLayout;
    footerLayout->setAlignment(Qt::AlignCenter); //can this be done in CSS???
    footerLayout->addWidget(nsfLogo);
    footerLayout->addWidget(nsfText);
    //footerLayout->addWidget(simLogo);

    footer->setLayout(footerLayout);

    largeLayout->addWidget(footer);
}

void MainWindow::createInputPanel() {
    inputLayout = new QVBoxLayout;


    //
    // Create a section line + title + add
    // styleSheet

    QFrame *line4 = new QFrame();
    line4->setObjectName(QString::fromUtf8("line"));
    line4->setGeometry(QRect(320, 150, 118, 3));
    line4->setFrameShape(QFrame::HLine);
    line4->setFrameShadow(QFrame::Sunken);

    QLabel *inTitle = new QLabel(); //styleSheet
    inTitle->setText(tr("Input Motion")); //styleSheet
    inTitle->setObjectName(QString::fromUtf8("inTitle")); //styleSheet

    inputLayout->addWidget(inTitle);
    inputLayout->addWidget(line4);


    //
    // create the frame for the input motion selection
    //

    QFrame *inputMotion = new QFrame(); //styleSheet
    inputMotion->setObjectName(QString::fromUtf8("inputMotion")); //styleSheet
    QHBoxLayout *inputMotionLayout = new QHBoxLayout();
    QLabel *sectionTitle = new QLabel();
    sectionTitle->setText(tr("Input Motion Title"));
    sectionTitle->setObjectName(QString::fromUtf8("sectionTitle")); //styleSheet
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(tr("Input Motion"));

    inMotion = new QComboBox();
    inputMotionLayout->addWidget(entryLabel);
    inputMotionLayout->addWidget(inMotion);
    addMotion = new QPushButton("Add");
    inputMotionLayout->addWidget(addMotion);
    inputMotion->setLayout(inputMotionLayout);
    // inputMotion->setFrameStyle(QFrame::Raised);
    inputMotion->setLineWidth(1);
    inputMotion->setFrameShape(QFrame::Box);

    inputLayout->addWidget(inputMotion);

    //
    // Create a section line
    // styleSheet

    QFrame *line = new QFrame();
    line->setObjectName(QString::fromUtf8("line"));
    line->setGeometry(QRect(320, 150, 118, 3));
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);


    //
    // Create a section line 2
    // styleSheet

    QFrame *line2 = new QFrame();
    line2->setObjectName(QString::fromUtf8("line2"));
    line2->setGeometry(QRect(320, 150, 118, 3));
    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);
    QLabel *propertiesTitle = new QLabel();
    propertiesTitle->setText(tr("Building Properties"));
    propertiesTitle->setObjectName(QString::fromUtf8("propertiesTitle"));
    inputLayout->addWidget(propertiesTitle);
    inputLayout->addWidget(line2);


    //
    // create to hold major model inputs
    //

    QFrame *mainProperties = new QFrame();
    mainProperties->setObjectName(QString::fromUtf8("mainProperties")); //styleSheet
    QVBoxLayout *mainPropertiesLayout = new QVBoxLayout();
    inFloors = createTextEntry(tr("number Floors"), mainPropertiesLayout);
    inWeight = createTextEntry(tr("building Weight"), mainPropertiesLayout);
    inHeight = createTextEntry(tr("building Height"), mainPropertiesLayout);
    inK = createTextEntry(tr("story Stiffness"), mainPropertiesLayout);
    inDamping = createTextEntry(tr("damping Ratio"), mainPropertiesLayout);
    inGravity =  createTextEntry(tr("gravity"), mainPropertiesLayout);
    pDeltaBox = new QCheckBox(tr("Include PDelta"), 0);
    pDeltaBox->setCheckState(Qt::Checked);

    mainPropertiesLayout->addWidget(pDeltaBox);
    mainProperties->setLayout(mainPropertiesLayout);
    // mainProperties->setFrameStyle(QFrame::Raised);
    mainProperties->setLineWidth(1);
    mainProperties->setFrameShape(QFrame::Box);
    inputLayout->addWidget(mainProperties);


    //
    // create frames to hold story and floor selections
    //


    floorMassFrame = new QFrame();
    floorMassFrame->setObjectName(QString::fromUtf8("floorMassFrame")); //styleSheet
    QVBoxLayout *floorMassFrameLayout = new QVBoxLayout();
    inFloorWeight = createTextEntry(tr("Floor Weight"), floorMassFrameLayout);
    floorMassFrame->setLayout(floorMassFrameLayout);
    floorMassFrame->setLineWidth(1);
    floorMassFrame->setFrameShape(QFrame::Box);
    inputLayout->addWidget(floorMassFrame);
    floorMassFrame->setVisible(false);

    storyPropertiesFrame = new QFrame();
    storyPropertiesFrame->setObjectName(QString::fromUtf8("storyPropertiesFrame"));
    QVBoxLayout *storyPropertiesFrameLayout = new QVBoxLayout();
    inStoryHeight = createTextEntry(tr("Story Height"), storyPropertiesFrameLayout);
    inStoryK = createTextEntry(tr("Stiffness"), storyPropertiesFrameLayout);
    inStoryFy = createTextEntry(tr("Yield Strength"), storyPropertiesFrameLayout);
    inStoryB = createTextEntry(tr("Hardening Ratio"), storyPropertiesFrameLayout);
    storyPropertiesFrame->setLayout(storyPropertiesFrameLayout);
    storyPropertiesFrame->setLineWidth(1);
    storyPropertiesFrame->setFrameShape(QFrame::Box);
    inputLayout->addWidget(storyPropertiesFrame);
    storyPropertiesFrame->setVisible(false);

    spreadSheetFrame = new QFrame();
    QVBoxLayout *spreadsheetFrameLayout = new QVBoxLayout();
    headings << tr("Weight") << tr("Heighht") << tr("K") << tr("Fy") << tr("b") << tr("zeta");
    dataTypes << SIMPLESPREADSHEET_QDouble;
    dataTypes << SIMPLESPREADSHEET_QDouble;
    dataTypes << SIMPLESPREADSHEET_QDouble;
    dataTypes << SIMPLESPREADSHEET_QDouble;
    dataTypes << SIMPLESPREADSHEET_QDouble;
    dataTypes << SIMPLESPREADSHEET_QDouble;

    theSpreadsheet = new SimpleSpreadsheetWidget(6, 4, headings, dataTypes,0);
    //theSpreadsheet = new QTableWidget();
    spreadsheetFrameLayout->addWidget(theSpreadsheet, 1.0);

    spreadSheetFrame->setLayout(spreadsheetFrameLayout);
    spreadSheetFrame->setLineWidth(1);
    spreadSheetFrame->setFrameShape(QFrame::Box);
    //theSpreadsheet->setS

    inputLayout->addWidget(spreadSheetFrame);
    spreadSheetFrame->setVisible(false);

    inputLayout->addStretch();

    //
    // finally create a frame to hold push buttons
    //

    QFrame *pushButtons = new QFrame();
    pushButtons->setObjectName(QString::fromUtf8("pushButtons")); //styleSheet
    QHBoxLayout *pushButtonsLayout = new QHBoxLayout();
    runButton = new QPushButton("Run");
    pushButtonsLayout->addWidget(runButton);
    stopButton = new QPushButton("Stop");
    pushButtonsLayout->addWidget(stopButton);
    exitButton = new QPushButton("Exit");
    pushButtonsLayout->addWidget(exitButton);
    pushButtons->setLayout(pushButtonsLayout);
    // mainProperties->setFrameStyle(QFrame::Raised);
    pushButtons->setLineWidth(1);
    pushButtons->setFrameShape(QFrame::Box);


    inputLayout->addWidget(pushButtons);

    mainLayout->addLayout(inputLayout);
    //largeLayout->addLayout(mainLayout);


    //
    // set validators for QlineEdits
    //

    inFloors->setValidator(new QIntValidator);
    inWeight->setValidator(new QDoubleValidator);
    inHeight->setValidator(new QDoubleValidator);
    inK->setValidator(new QDoubleValidator);
    inDamping->setValidator(new QDoubleValidator);
    inGravity->setValidator(new QDoubleValidator);
    inFloorWeight->setValidator(new QDoubleValidator);
    inStoryB->setValidator(new QDoubleValidator);
    inStoryFy->setValidator(new QDoubleValidator);
    inStoryHeight->setValidator(new QDoubleValidator);
    inStoryK->setValidator(new QDoubleValidator);

    //
    // connect signals & slots
    //
    connect(pDeltaBox, SIGNAL(stateChanged(int)), this, SLOT(on_includePDeltaChanged(int)));
    connect(addMotion,SIGNAL(clicked()), this, SLOT(on_addMotion_clicked()));
    connect(inFloors,SIGNAL(editingFinished()), this, SLOT(on_inFloors_editingFinished()));
    connect(inWeight,SIGNAL(editingFinished()), this, SLOT(on_inWeight_editingFinished()));
    connect(inHeight,SIGNAL(editingFinished()), this, SLOT(on_inHeight_editingFinished()));
    connect(inK,SIGNAL(editingFinished()), this, SLOT(on_inK_editingFinished()));
    connect(inDamping,SIGNAL(editingFinished()), this, SLOT(on_inDamping_editingFinished()));
    connect(inFloorWeight,SIGNAL(editingFinished()), this, SLOT(on_inFloorWeight_editingFinished()));
    connect(inStoryB, SIGNAL(editingFinished()),this,SLOT(on_inStoryB_editingFinished()));
    connect(inStoryFy,SIGNAL(editingFinished()),this, SLOT(on_inStoryFy_editingFinished()));
    connect(inStoryHeight, SIGNAL(editingFinished()), this, SLOT(on_inStoryHeight_editingFinished()));
    connect(inStoryK, SIGNAL(editingFinished()), this, SLOT(on_inStoryK_editingFinished()));

    connect(inMotion, SIGNAL(currentIndexChanged(QString)), this, SLOT(on_inMotionSelection_currentTextChanged(QString)));

    connect(theSpreadsheet, SIGNAL(cellClicked(int,int)), this, SLOT(on_theSpreadsheet_cellClicked(int,int)));
    connect(theSpreadsheet, SIGNAL(cellEntered(int,int)), this,SLOT(on_theSpreadsheet_cellClicked(int,int)));
    connect(theSpreadsheet, SIGNAL(cellChanged(int,int)), this,SLOT(on_theSpreadsheet_cellChanged(int,int)));

    connect(runButton, SIGNAL(clicked()), this, SLOT(on_runButton_clicked()));
    connect(stopButton, SIGNAL(clicked()), this, SLOT(on_stopButton_clicked()));
    connect(exitButton, SIGNAL(clicked()), this, SLOT(on_exitButton_clicked()));
}

void MainWindow::createOutputPanel() {

    // 1) Basic Outputs, e.g. Disp, Periods
    // 2) MyGlWidget
    // 3) QCustomPlotWidget
    // 4) CurrentTime

    outputLayout = new QVBoxLayout;

    //
    // Create a section line + title + add
    // styleSheet

    QFrame *line3 = new QFrame();
    line3->setObjectName(QString::fromUtf8("line"));
    line3->setGeometry(QRect(320, 150, 118, 3));
    line3->setFrameShape(QFrame::HLine);
    line3->setFrameShadow(QFrame::Sunken);

    QLabel *outTitle = new QLabel(); //styleSheet
    outTitle->setText(tr("Output")); //styleSheet
    outTitle->setObjectName(QString::fromUtf8("outTitle")); //styleSheet

    outputLayout->addWidget(outTitle);
    outputLayout->addWidget(line3);


    // frame for basic outputs,
    QFrame *outputMaxFrame = new QFrame();

    // frame for max disp / current period
    QFrame *firstOutput = new QFrame(); //styleSheet
    firstOutput->setObjectName(QString::fromUtf8("firstOutput"));
    QVBoxLayout *firstOutputLayout = new QVBoxLayout();
    maxDispLabel = createLabelEntry(tr("Max Disp"), firstOutputLayout); //styleSheet
    currentPeriod= createLabelEntry(tr("Fundamental Period"),firstOutputLayout); //styleSheet
    firstOutput->setLayout(firstOutputLayout);
    firstOutput->setLineWidth(1);
    firstOutput->setFrameShape(QFrame::Box);
    outputLayout->addWidget(firstOutput);


    QVBoxLayout *outputMaxLayout = new QVBoxLayout();
    QLabel *vizTitle = new QLabel(); //styleSheet
    vizTitle->setText(tr("Visualization Section Title")); //styleSheet
    vizTitle->setObjectName(QString::fromUtf8("vizTitle")); //styleSheet
    //maxDispLabel = createLabelEntry(tr("Max Disp"), firstOutputLayout); //styleSheet
    //currentPeriod= createLabelEntry(tr("Fundamental Period"),firstOutputLayout); //styleSheet
    outputMaxFrame->setLayout(outputMaxLayout);
    outputMaxFrame->setLineWidth(1);
    outputMaxFrame->setFrameShape(QFrame::Box);
    outputMaxFrame->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    outputMaxFrame->setLayout(firstOutputLayout); //this does not set properly???????
    //outputLayout->addWidget(outputMaxFrame);

    //
    // Create a section line
    // styleSheet

    QFrame *line = new QFrame();
    line->setObjectName(QString::fromUtf8("line"));
    line->setGeometry(QRect(320, 150, 118, 3));
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);

    //
    // Add title and line
    // styleSheet

    //outputLayout->addWidget(vizTitle);
    //outputLayout->addWidget(line);

    // GL Widget
    myGL = new MyGlWidget();
    myGL->setMinimumHeight(400);
    myGL->setMinimumWidth(250);
    myGL->setModel(this);
    outputLayout->addWidget(myGL);

    // input acceleration plot
    earthquakePlot=new QCustomPlot();
    earthquakePlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);
    earthquakePlot->setMinimumHeight(100);
    earthquakePlot->setMaximumHeight(100);
    outputLayout->addWidget(earthquakePlot);

    // slider for manual movement
    slider=new QSlider(Qt::Horizontal);
    outputLayout->addWidget(slider);

    // output frame to show current time
    QFrame *outputDataFrame = new QFrame();
    outputDataFrame->setObjectName(QString::fromUtf8("outputDataFrame"));
    QVBoxLayout *outputDataLayout = new QVBoxLayout();
    currentTime = createLabelEntry(tr("Current Time"), outputDataLayout);
    currentDisp = createLabelEntry(tr("Current Roof Disp"), outputDataLayout);
    outputDataFrame->setLayout(outputDataLayout);
    outputDataFrame->setLineWidth(1);
    outputDataFrame->setFrameShape(QFrame::Box);
    outputDataFrame->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    outputLayout->addWidget(outputDataFrame);

    // add layout to mainLayout and to largeLayout
    mainLayout->addLayout(outputLayout);
    largeLayout->addLayout(mainLayout); //styleSheet

    // signal and slot connects for slider
    connect(slider, SIGNAL(sliderPressed()),  this, SLOT(on_slider_sliderPressed()));
    connect(slider, SIGNAL(sliderReleased()), this, SLOT(on_slider_sliderReleased()));
    connect(slider, SIGNAL(valueChanged(int)),this, SLOT(on_slider_valueChanged(int)));
}

