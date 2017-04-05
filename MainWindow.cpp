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

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
Domain theDomain;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    numFloors(0), period(0), buildingW(0),storyK(0),
    weights(0), k(0), fy(0), b(0), floorHeights(0), storyHeights(0),
    dampingRatio(0.02), g(386.4), dt(0), gMotion(0),
    needAnalysis(true), eqData(0), dispResponses(0), maxDisp(0),
    movingSlider(false), fMinSelected(-1),fMaxSelected(-1), sMinSelected(-1),sMaxSelected(-1),
    time(1560),values(1560), graph(0), groupTracer(0)
{
    ui->setupUi(this);
    //ui->numFloors->setValidator( new QIntValidator);

    ui->inDamping->setValidator(new QDoubleValidator);
    ui->inFloors->setValidator(new QIntValidator);
    ui->inHeight->setValidator(new QDoubleValidator);
    ui->inK->setValidator(new QDoubleValidator);
    ui->inPeriod->setValidator(new QDoubleValidator);
    ui->myGL->setModel(this);
    ui->inFloorWeight->setValidator(new QDoubleValidator);
    ui->inGravity->setValidator(new QDoubleValidator);
    ui->inStoryB->setValidator(new QDoubleValidator);
    ui->inStoryFy->setValidator(new QDoubleValidator);
    ui->inStoryHeight->setValidator(new QDoubleValidator);
    ui->inStoryK->setValidator(new QDoubleValidator);

    // create elCentro EarthquakeRecord and make current
    QStringList elCentrolist = elCentroTextData.split(QRegExp("[\r\n\t ]+"), QString::SkipEmptyParts);
    Vector *elCentroData = new Vector(elCentrolist.size()+1);
    (*elCentroData)(0) = 0;
    time[0]=0.;
    values[0]=0.;
    double maxValue = 0;
    for (int i = 0; i < elCentrolist.size(); ++i) {
        double value = elCentrolist.at(i).toDouble();
        (*elCentroData)(i+1) = value;
        time[i+1]=i*0.02;
        values[i+1]=value;
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
    QString blankString("BLANK");
    EarthquakeRecord *blank = new EarthquakeRecord(blankString, 100, 0.02, blankData);
    records.insert(std::make_pair(blankString, blank));

    this->setBasicModel(4, 0.4);

    ui->inFloorWeight->setDisabled(true);

    ui->inStoryHeight->setDisabled(true);
    ui->inStoryK->setDisabled(true);
    ui->inStoryB->setDisabled(true);
    ui->inStoryFy->setDisabled(true);

    ui->inHazard->addItem(QString("Earthquake"));



    ui->inMotionSelection->addItem(elCentroString);
    ui->inMotionSelection->addItem(QString("BLANK"));

    ui->slider->setRange(0, numSteps);
    ui->slider->setSliderPosition(0);
    //ui->slider->setMaximum(numSteps);

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
}

void MainWindow::draw(MyGlWidget *theGL)
{
    if (needAnalysis == true) {
        doAnalysis();
    }

    for (int i=0; i<numFloors; i++) {
 if (i >= sMinSelected && i <= sMaxSelected)
        theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                        dispResponses[i+1][currentStep],floorHeights[i+1], 2, 1, 0, 0);
    else
        theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                        dispResponses[i+1][currentStep],floorHeights[i+1], 2, 0, 0, 0);
    }

    for (int i=0; i<=numFloors; i++) {
       if (i >= fMinSelected && i <= fMaxSelected)
        theGL->drawNode(i, dispResponses[i][currentStep],floorHeights[i], 10, 1, 0, 0);
       else
        theGL->drawNode(i, dispResponses[i][currentStep],floorHeights[i], 10, 0, 0, 1);
    }

    // display range of displacement
    static char maxDispString[30];
    snprintf(maxDispString, 50, "%.3e", maxDisp);
    theGL->drawLine(0, -maxDisp, 0.0, maxDisp, 0.0, 1.0, 0., 0., 0.);
    theGL->drawText(0, -maxDisp, buildingH/100., maxDispString,0,0,0);
    //theGL->drawText(0, maxDisp, buildingH/100., maxDispString,0,0,0);

    // display current time

    ui->currentTime->setText(QString().setNum(currentStep*dt,'f',2));

    // update red dot on earthquake plot
    groupTracer->setGraph(0);
    groupTracer->setGraph(graph);
    groupTracer->setGraphKey(currentStep*dt);
    groupTracer->updatePosition();
    ui->earthquakePlot->replot();
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

        dispResponses = new double *[numF+1];
        for (int i=0; i<numF+1; i++) {
            dispResponses[i] = new double[numSteps+1]; // +1 as doing 0 at start
        }
    }

    // set values
    double floorW = W/(numF);

    for (int i=0; i<numF; i++) {
      weights[i] = floorW;
      k[i] = K;
      fy[i] = 1.0e100;
      b[i] = 0.;
      floorHeights[i] = i;
      storyHeights[i] = 1;
    }
    floorHeights[numF] = numF;

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
    ui->inGravity->setText(QString::number(g));
    needAnalysis = true;
}

void MainWindow::on_inFloors_editingFinished()
{
    QString textFloors =  ui->inFloors->text();
    int numFloorsText = textFloors.toInt();
    if (numFloorsText != numFloors) {
        this->setBasicModel(numFloorsText, period);
    }
  //  ui->inWeight->setFocus();
    this->reset();
}

void MainWindow::on_inWeight_editingFinished()
{
    QString textW =  ui->inWeight->text();
    double textToDoubleW = textW.toDouble();
    if (textToDoubleW != buildingW) {
        // set values
        double floorW = textToDoubleW/(numFloors);
        for (int i=0; i<numFloors; i++) {
            weights[i] = floorW;
        }
    }
    // ui->inHeight->setFocus();
    this->reset();
}

void MainWindow::on_inHeight_editingFinished()
{
    QString textH =  ui->inHeight->text();
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
  //   ui->inK->setFocus();
    this->reset();
}


void MainWindow::on_inK_editingFinished()
{
    QString text =  ui->inK->text();
    double textToDouble = text.toDouble();
    for (int i=0; i<=numFloors; i++)
        k[i] = textToDouble;
  //  ui->inDamping->setFocus();
    this->reset();
}

void MainWindow::on_inDamping_editingFinished()
{
    QString text =  ui->inDamping->text();
    double textToDouble = text.toDouble();
    dampingRatio = textToDouble;
//    ui->inGravity->setFocus();
    this->reset();
}

void MainWindow::on_inFloorWeight_editingFinished()
{
    QString text =  ui->inFloorWeight->text();
    double textToDouble = text.toDouble();
    for (int i=fMinSelected; i<=fMaxSelected; i++)
        weights[i] = textToDouble;

    this->reset();
}




void MainWindow::on_inGravity_editingFinished()
{
    QString text =  ui->inGravity->text();
    double textToDouble = text.toDouble();
    g = textToDouble;

    this->reset();
}

void MainWindow::on_inStoryHeight_editingFinished()
{
    QString text =  ui->inStoryHeight->text();
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
    delete [] floorHeights;
    floorHeights = newFloorHeights;

    // move focus, update graphic and set analysis flag
    ui->inStoryK->setFocus();
    this->reset();
}

void MainWindow::on_inStoryK_editingFinished()
{
    QString text =  ui->inStoryK->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        k[i] = textToDouble;

    ui->inStoryFy->setFocus();
    this->reset();
}

void MainWindow::on_inStoryFy_editingFinished()
{
    QString text =  ui->inStoryFy->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        fy[i] = textToDouble;

    ui->inStoryB->setFocus();
    this->reset();
}

void MainWindow::on_inStoryB_editingFinished()
{
    QString text =  ui->inStoryB->text();
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        b[i] = textToDouble;

    ui->inStoryHeight->setFocus();
    this->reset();
   // needAnalysis = true;
   // ui->myGL->update();
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
                theMass(0,0) = weights[i-1]/g;
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
        EquiSolnAlgo      *theSolnAlgo = new NewtonRaphson();
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

        //
        //analyze & get results
        //
        maxDisp = 0;
        for (int i=0; i<=numSteps; i++) { // <= due to adding 0 at start
            theAnalysis.analyze(1, dt);
            for (int j=0; j<numFloors+1; j++) {
                double nodeDisp = theNodes[j]->getDisp()(0);
                dispResponses[j][i] = nodeDisp;
                if (fabs(nodeDisp) > maxDisp)
                        maxDisp = fabs(nodeDisp);
            }
        }

        // clean up memory
        delete [] theNodes;

        // reset values, i.e. slider position, current displayed step, and display properties
        needAnalysis = false;
        currentStep = 0;
        groupTracer->setGraphKey(0);
        ui->slider->setSliderPosition(0);
        ui->myGL->update();
    }
}

void MainWindow::reset() {
    needAnalysis = true;
    ui->myGL->update();

    // update the properties table
    ui->tableWidget->clear();
    ui->tableWidget->setColumnCount(5);
    ui->tableWidget->setRowCount(numFloors);
    ui->tableWidget->horizontalHeader()->setStretchLastSection(true);// horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
   // ui->tableWidget->setFixedWidth(344);
    updatingPropertiesTable = true;
    ui->tableWidget->setHorizontalHeaderLabels(QString(" Weight ; Height ;    K    ;    Fy    ;    b    ").split(";"));
    for (int i=0; i<numFloors; i++) {
        ui->tableWidget->setItem(i,0,new QTableWidgetItem(QString().setNum(weights[i])));
        ui->tableWidget->setItem(i,1,new QTableWidgetItem(QString().setNum(storyHeights[i])));
        ui->tableWidget->setItem(i,2,new QTableWidgetItem(QString().setNum(k[i])));
        ui->tableWidget->setItem(i,3,new QTableWidgetItem(QString().setNum(fy[i])));
        ui->tableWidget->setItem(i,4,new QTableWidgetItem(QString().setNum(b[i])));
    }
    ui->tableWidget->resizeRowsToContents();
    ui->tableWidget->resizeColumnsToContents();

    updatingPropertiesTable = false;
}

void MainWindow::on_tableWidget_cellChanged(int row, int column)
{
    if (updatingPropertiesTable == false) {
    QString text = ui->tableWidget->item(row,column)->text();
    bool ok;
     double textToDouble = text.toDouble(&ok);
     if (column == 0) {
         weights[row] = textToDouble;
     } else if (column  == 1) {
         storyHeights[row] = textToDouble;
         for (int i=row; i<numFloors; i++)
             floorHeights[i+1] = floorHeights[i]+storyHeights[i];
         buildingH = floorHeights[numFloors];
    } else if (column == 2) {
        k[row] = textToDouble;
     } else if (column == 3) {
         fy[row] = textToDouble;
     } else
         b[row] = textToDouble;
    }
    needAnalysis = true;
    ui->myGL->update();
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
    do { //while (currentStep < numSteps && stopRun == false){
        ui->slider->setSliderPosition(currentStep);
        ui->myGL->repaint();
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
            ui->myGL->update();
        }
        currentStep = ui->slider->value();

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
    qDebug() << numFloors;
    for (int i=0; i<=numFloors; i++) {
        if (floorHeights[i] >= yMin && floorHeights[i] <= yMax) {
            qDebug() << i << " " << floorHeights[i] << " " << yMin << " " <<  yMax;
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
    qDebug() << "sMinSelected: " << sMinSelected << " sMaxSelected: " << sMaxSelected;
    //
    // based on min, max nodes enable/disable lineEdits & set text
    //
    if (fMinSelected == -1) {
        // blank out line edits
        QString blank("");
        //ui->inFloorHeight->setText(blank);
        // disable line edits
        ui->inFloorWeight->setDisabled(true);

    } else {

        ui->inFloorWeight->setDisabled(false);

    }
    if (sMaxSelected == -1) { // (fMaxSelected < fMinSelected) {
        // blank out line edits
        QString blank("");
        ui->inStoryHeight->setText(blank);
        ui->inStoryK->setText(blank);
        ui->inStoryB->setText(blank);
        ui->inStoryFy->setText(blank);

        ui->inStoryHeight->setDisabled(true);
        ui->inStoryK->setDisabled(true);
        ui->inStoryB->setDisabled(true);
        ui->inStoryFy->setDisabled(true);

    } else {
        // enable line edits

        ui->inStoryHeight->setDisabled(false);
        ui->inStoryK->setDisabled(false);
        ui->inStoryB->setDisabled(false);
        ui->inStoryFy->setDisabled(false);

        // set text line edits
    }
    ui->myGL->repaint();
}




void MainWindow::on_tableWidget_cellClicked(int row, int column)
{
    qDebug() << row << " " << column;
    fMinSelected = row+1; fMaxSelected = row+1;
    sMinSelected = row; sMaxSelected = row;
    ui->myGL->repaint();
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
        values.resize(numSteps);
        time.resize(numSteps);
        for (int i = 0; i < numSteps; ++i) {
            double value = (*eqData)[i];
            time[i]=i*dt;
            values[i]=value;
            if (fabs(value) > maxValue)
                maxValue = fabs(value);
        }

        // reset earthquake plot
        ui->earthquakePlot->clearGraphs();
        graph = ui->earthquakePlot->addGraph();
        ui->earthquakePlot->graph(0)->setData(time, values);
        ui->earthquakePlot->xAxis->setRange(0, numSteps*dt);
        ui->earthquakePlot->yAxis->setRange(-maxValue, maxValue);
        ui->earthquakePlot->axisRect()->setAutoMargins(QCP::msNone);
        ui->earthquakePlot->axisRect()->setMargins(QMargins(0,0,0,0));
        if (groupTracer != 0)
            delete groupTracer;
        groupTracer = new QCPItemTracer(ui->earthquakePlot);
        groupTracer->setGraph(graph);
        groupTracer->setGraphKey(0);
        groupTracer->setInterpolating(true);
        groupTracer->setStyle(QCPItemTracer::tsCircle);
        groupTracer->setPen(QPen(Qt::red));
        groupTracer->setBrush(Qt::red);
        groupTracer->setSize(7);

        // reset slider range
        ui->slider->setRange(0, numSteps);


        this->reset();
        qDebug() << numSteps << " " << dt;
    }
}




void MainWindow::on_pushButton_2_released()
{
  QApplication::quit();
}
