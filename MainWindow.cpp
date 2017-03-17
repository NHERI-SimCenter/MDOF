#include "MainWindow.h"
#include "ui_MainWindow.h"

#include <QDebug>
#include <QSlider>


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
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <DirectIntegrationAnalysis.h>
#include <AnalysisModel.h>
#include "elCentro.AT2"
#include <Vector.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
Domain theDomain;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    numFloors(0), period(0), buildingW(0),storyK(0),
    masses(0), k(0), fy(0), b(0), heights(0),
    dampingRatio(0.02), dt(0), gMotion(0),
    needAnalysis(true), elCentroData(0), dispResponses(0), maxDisp(0),
    movingSlider(false), fMinSelected(-1),fMaxSelected(-1)
{
    ui->setupUi(this);
    //ui->numFloors->setValidator( new QIntValidator);

    ui->inDamping->setValidator(new QDoubleValidator);
    ui->inFloors->setValidator(new QIntValidator);
    ui->inHeight->setValidator(new QDoubleValidator);
    ui->inK->setValidator(new QDoubleValidator);
    ui->inPeriod->setValidator(new QDoubleValidator);
    ui->myGL->setModel(this);
    this->setBasicModel(4, 0.4);

    ui->inFloorWeight->setDisabled(true);
    ui->inFloorHeight->setDisabled(true);

    ui->inStoryHeight->setDisabled(true);
    ui->inStoryK->setDisabled(true);
    ui->inStoryB->setDisabled(true);
    ui->inStoryFy->setDisabled(true);

    QStringList elCentrolist = elCentroTextData.split(QRegExp("[\r\n\t ]+"), QString::SkipEmptyParts);
    elCentroData = new Vector(elCentrolist.size()+1);
    (*elCentroData)(0) = 0;
    for (int i = 0; i < elCentrolist.size(); ++i) {
        (*elCentroData)(i+1) = elCentrolist.at(i).toDouble();
    }
    dt = 0.02;
    numSteps = 1560;
    ui->slider->setRange(0, numSteps);
    ui->slider->setSliderPosition(0);
    //ui->slider->setMaximum(numSteps);
}

MainWindow::~MainWindow()
{
    delete ui;

    if (masses != 0)
        delete [] masses;
    if (k != 0)
        delete [] k;
    if (fy != 0)
        delete [] fy;
    if (b != 0)
        delete [] b;
    if (gMotion != 0)
        delete [] gMotion;
}

void MainWindow::draw(MyGlWidget *theGL)
{
    if (needAnalysis == true) {
        doAnalysis();
    }

    for (int i=0; i<numFloors; i++) {
 if (i >= fMinSelected && i+1 <= fMaxSelected)
        theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],heights[i],
                        dispResponses[i+1][currentStep],heights[i+1], 2, 1, 0, 0);
    else
        theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],heights[i],
                        dispResponses[i+1][currentStep],heights[i+1], 2, 0, 0, 0);
    }

    for (int i=0; i<=numFloors; i++) {
       if (i >= fMinSelected && i <= fMaxSelected)
        theGL->drawNode(i, dispResponses[i][currentStep],heights[i], 10, 1, 0, 0);
       else
        theGL->drawNode(i, dispResponses[i][currentStep],heights[i], 10, 0, 0, 1);
    }
    ui->currentTime->setText(QString().setNum(currentStep*dt,'f',2));

    if (currentStep < numSteps)
       currentStep++;
}



void MainWindow::updatePeriod()
{
    period = 2.0;
}

void MainWindow::setBasicModel(int numF, double period)
{
    storyK = 1.0;
    buildingW = 1.0; // FMK calculate based on period

    this->setBasicModel(numF, buildingW, storyK);
}

void MainWindow::setBasicModel(int numF, double W, double K)
{
    if (numFloors != numF) {
        // if invalid numFloor, return
        if (numF <= 0)
            return;

        // resize arrays
        if (masses != 0)
            delete [] masses;
        if (k != 0)
            delete [] k;
        if (fy != 0)
            delete [] fy;
        if (b != 0)
            delete [] b;
        if (dispResponses != 0) {
               for (int j=0; j<numFloors+1; j++)
                   delete [] dispResponses[j];
               delete [] dispResponses;
        }

        masses = new double[numF];
        k = new double[numF];
        fy = new double[numF];
        b = new double[numF];
        heights = new double[numF+1];

        dispResponses = new double *[numF+1];
        for (int i=0; i<numF+1; i++) {
            dispResponses[i] = new double[numSteps];
        }
    }

    // set values
    double grav = 386.4;
    double floorM = W/(numF*grav);

    for (int i=0; i<numF; i++) {
        masses[i] = floorM;
        k[i] = K;
        fy[i] = 1.0e100;
        b[i] = 0.;
        heights[i] = i;
    }
    heights[numF] = numF;

    buildingW = W;
    storyK = K;
    numFloors = numF;
    buildingH = numF;

    this->updatePeriod();

    // update text boxes
    ui->inPeriod->setText(QString::number(period));
    ui->inWeight->setText(QString::number(buildingW));
    ui->inK->setText(QString::number(storyK));
    ui->inHeight->setText(QString::number(numF));

    needAnalysis = true;
}

void MainWindow::on_inFloors_editingFinished()
{
    QString textFloors =  ui->inFloors->text();
    int numFloorsText = textFloors.toInt();
    if (numFloorsText != numFloors) {
        this->setBasicModel(numFloorsText, period);
    }

    needAnalysis = true;
    ui->myGL->update();
}

void MainWindow::on_inWeight_editingFinished()
{
    QString textW =  ui->inWeight->text();
    double textToDoubleW = textW.toDouble();
    if (textToDoubleW != buildingW) {
        // set values
        double grav = 386.4;
        double floorM = textToDoubleW/(numFloors*grav);
        for (int i=0; i<numFloors; i++) {
            masses[i] = floorM;
        }
    }
    needAnalysis = true;
    ui->myGL->update();
}

void MainWindow::on_inHeight_editingFinished()
{
    QString textH =  ui->inHeight->text();
    double textToDoubleH = textH.toDouble();
    if (textToDoubleH != buildingH) {
        // set values
        buildingH = textToDoubleH;
        double deltaH = buildingH/numFloors;
        // qDebug() << deltaH << " " << buildingH << " " << numFloors;
        heights[0] = 0;
        for (int i=1; i<=numFloors; i++) {
            heights[i] = deltaH + heights[i-1];
        }
    }
    needAnalysis = true;
    ui->myGL->update();
}


void MainWindow::doAnalysis()
{
    if (needAnalysis == true) {

        // clear existing model
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
                theMass(0,0) = masses[i];
                theNode->setMass(theMass);
            }
        }


        // create the vectors for the element orientation
        Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
        Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

        for (int i=0; i<numFloors; i++) {
            UniaxialMaterial *theMat = new Steel01(i+1,fy[i],k[i],b[i]);
            ZeroLength *theEle = new ZeroLength(i+1, 1, i+1, i+2,
                                                x, y, *theMat, 0);
            theDomain.addElement(theEle);
        }

        //
        // create load pattern and add loads
        //

        PathSeries *theSeries = new PathSeries(1, *elCentroData, dt, 386.4);
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
        EquiSolnAlgo      *theSolnAlgo = new NewtonRaphson();
        TransientIntegrator  *theIntegrator = new Newmark(0.5, 0.25);
        //ConstraintHandler *theHandler = new TransformationConstraintHandler();
        ConstraintHandler *theHandler = new PlainHandler();
        RCM               *theRCM = new RCM();
        DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);
        BandGenLinSolver  *theSolver = new BandGenLinLapackSolver();
        LinearSOE         *theSOE = new BandGenLinSOE(*theSolver);

        DirectIntegrationAnalysis    theAnalysis(theDomain,
                                                 *theHandler,
                                                 *theNumberer,
                                                 *theModel,
                                                 *theSolnAlgo,
                                                 *theSOE,
                                                 *theIntegrator);
        theSolnAlgo->setConvergenceTest(theTest);

        //
        //analyze & get results
        //
        maxDisp = 0;
        for (int i=0; i<numSteps; i++) {
            theAnalysis.analyze(1, dt);
            for (int j=0; j<numFloors+1; j++) {
                double nodeDisp = theNodes[j]->getDisp()(0);
                dispResponses[j][i] = nodeDisp;
                if (fabs(nodeDisp) > maxDisp)
                        maxDisp = fabs(nodeDisp);
            }
        }

        // reset values, i.e. slider position, current displayed step, and display properties
        needAnalysis = false;
        currentStep = 0;
        ui->slider->setSliderPosition(0);
        ui->myGL->update();
    }
}


void MainWindow::on_stopButton_clicked()
{
    stopRun = true;
}

void MainWindow::on_runButton_clicked()
{
    stopRun = false;
    if (needAnalysis == true) {
        this->doAnalysis();

    }

    currentStep = 0;
    while (currentStep < numSteps && stopRun == false){
        ui->slider->setSliderPosition(currentStep);
        ui->myGL->repaint();
        QCoreApplication::processEvents();
    }
}


void MainWindow::on_slider_valueChanged(int value)
{
    if (movingSlider == true) {
        stopRun = true;
        if (needAnalysis == true) {
            this->doAnalysis();
            ui->myGL->update();
        }
        currentStep = ui->slider->value();
        //qDebug() << currentStep;
        ui->myGL->repaint();
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
    float yMin = 0;
    float yMax = 0;
    if (y1 < y2) {
        yMin = y1;
        yMax = y2;
    } else {
        yMin = y2;
        yMax = y1;
    }
    fMaxSelected = -1;
    for (int i=0; i<numFloors+1; i++) {
        if (heights[i] < yMax)
            fMaxSelected = i;
    }
    fMinSelected = numFloors+2;
    for (int i=numFloors+1; i>=0; i--) {
        if (heights[i] > yMin)
            fMinSelected = i;
    }

    ui->myGL->repaint();
}
