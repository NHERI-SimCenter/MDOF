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
 OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
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
#include "sectiontitle.h"

#include <HeaderWidget.h>
#include <FooterWidget.h>


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
#include <ResponseWidget.h>

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
                int maxL=100,
                QString *unitText =0)
{
    QHBoxLayout *entryLayout = new QHBoxLayout();
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(text);

    QLineEdit *res = new QLineEdit();
    res->setMinimumWidth(minL);
    res->setMaximumWidth(maxL);
    res->setValidator(new QDoubleValidator);

    entryLayout->addWidget(entryLabel);
    entryLayout->addStretch();
    entryLayout->addWidget(res);

    if (unitText != 0) {
        QLabel *unitLabel = new QLabel();
        unitLabel->setText(*unitText);
        unitLabel->setMinimumWidth(40);
        unitLabel->setMaximumWidth(50);
        entryLayout->addWidget(unitLabel);

    }

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
                 int maxL=100,
                 QString *unitText = 0)
{
    QHBoxLayout *entryLayout = new QHBoxLayout();
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(text);

    QLabel *res = new QLabel();
    res->setMinimumWidth(minL);
    res->setMaximumWidth(maxL);
    res->setAlignment(Qt::AlignRight);

    entryLayout->addWidget(entryLabel);
    entryLayout->addStretch();
    entryLayout->addWidget(res);

    if (unitText != 0) {
        QLabel *unitLabel = new QLabel();
        unitLabel->setText(*unitText);
        unitLabel->setMinimumWidth(40);
        unitLabel->setMaximumWidth(100);
        entryLayout->addWidget(unitLabel);

    }
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
    includePDelta(true), needAnalysis(true), analysisFailed(false), motionData(0),
    dispResponses(0), storyForceResponses(0), storyDriftResponses(0), maxDisp(1),
    movingSlider(false), fMinSelected(-1),fMaxSelected(-1), sMinSelected(-1),sMaxSelected(-1),
    time(1561),excitationValues(1561), graph(0), groupTracer(0),floorSelected(-1),storySelected(-1)
{

    scaleFactor = 1.0;
    eigValues = new Vector;
    createActions();

    // create a main layout
    mainLayout = new QHBoxLayout();
    largeLayout = new QVBoxLayout();

    // create input and output panel layouts and each to main layout
    createHeaderBox();
    createInputPanel();
    createOutputPanel();
    createFooterBox();


    // create a widget in which to show everything //ALSO SET TO LARGE LAYOUT
    QWidget *widget = new QWidget();
    widget->setLayout(largeLayout);
    this->setCentralWidget(widget);

    QRect rec = QApplication::desktop()->screenGeometry();

    int height = 0.7*rec.height();
    int width = 0.7*rec.width();

    this->resize(width, height);

    //
    // create 2 blank motions & make elCentro current
    //

    QFile file(":/images/ElCentro.json");
    if(file.open(QFile::ReadOnly)) {
       QString jsonText = QLatin1String(file.readAll());
       QJsonDocument jsonDoc = QJsonDocument::fromJson(jsonText.toUtf8());
       QJsonObject jsonObj = jsonDoc.object();
       EarthquakeRecord *elCentro = new EarthquakeRecord();
       elCentro->inputFromJSON(jsonObj);

       QString recordString("ElCentro");
       records.insert(std::make_pair(QString("ElCentro"), elCentro));
       eqMotion->addItem(recordString);
    }

    QFile fileR(":/images/Rinaldi.json");
    if(fileR.open(QFile::ReadOnly)) {
       QString jsonText = QLatin1String(fileR.readAll());
       QJsonDocument jsonDoc = QJsonDocument::fromJson(jsonText.toUtf8());
       QJsonObject jsonObj = jsonDoc.object();
       EarthquakeRecord *rinaldi = new EarthquakeRecord();
       rinaldi->inputFromJSON(jsonObj);

       QString recordString("Northridge-Rinaldi");
       records.insert(std::make_pair(QString("Northridge-Rinaldi"), rinaldi));
       eqMotion->addItem(recordString);
    }

    // create a basic model with defaults
    this->setBasicModel(5, 5*100, 5*144, 31.54, .05, 386.4);

    // access a web page which will increment the usage count for this tool
    manager = new QNetworkAccessManager(this);

    connect(manager, SIGNAL(finished(QNetworkReply*)),
            this, SLOT(replyFinished(QNetworkReply*)));

    manager->get(QNetworkRequest(QUrl("http://opensees.berkeley.edu/OpenSees/developer/mdof/use.php")));
    manager->get(QNetworkRequest(QUrl("https://simcenter.designsafe-ci.org/multiple-degrees-freedom-analytics/")));
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
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 1, 0, 0);
        else if (i >= sMinSelected && i <= sMaxSelected)
            theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 1, 0, 0);
        else
            theGL->drawLine(i+1+numFloors, dispResponses[i][currentStep],floorHeights[i],
                            dispResponses[i+1][currentStep],floorHeights[i+1], 2, 0, 0, 0);
    }
    
    for (int i=0; i<=numFloors; i++) {
        if (i == floorSelected)
            theGL->drawPoint(i, dispResponses[i][currentStep],floorHeights[i], 10, 1, 0, 0);
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
        if (numF <= 0) {
            numF = 1;
        }

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
        if (storyForceResponses != 0) {
            for (int j=0; j<numFloors; j++)
                delete [] storyForceResponses[j];
            delete [] storyForceResponses;
        }
        if (storyDriftResponses != 0) {
            for (int j=0; j<numFloors; j++)
                delete [] storyDriftResponses[j];
            delete [] storyDriftResponses;
        }

        weights = new double[numF];
        k = new double[numF];
        fy = new double[numF];
        b = new double[numF];
        floorHeights = new double[numF+1];
        storyHeights = new double[numF];
        dampRatios = new double[numF];

        dispResponses = new double *[numF+1];
        storyForceResponses = new double *[numF];
        storyDriftResponses = new double *[numF];

        for (int i=0; i<numF+1; i++) {
            dispResponses[i] = new double[numSteps+1]; // +1 as doing 0 at start
            if (i<numF) {
                storyForceResponses[i] = new double[numSteps+1];
                storyDriftResponses[i] = new double[numSteps+1];
            }
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

    scaleFactorEQ->setText(QString::number(scaleFactor));
    inWeight->setText(QString::number(buildingW));
    inK->setText(QString::number(storyK));
    inFloors->setText(QString::number(numF));
    inHeight->setText(QString::number(buildingH));
    inDamping->setText(QString::number(zeta));

    periodHarmonicMotion = 0.5;
    magHarmonicMotion = 1.0;
    dtHarmonicMotion = 0.02;
    tFinalHarmonicMotion = 10.0;

    periodHarmonic->setText(QString::number(periodHarmonicMotion));
    magHarmonic->setText(QString::number(magHarmonicMotion));
    dtHarmonic->setText(QString::number(dtHarmonicMotion));
    tFinalHarmonic->setText(QString::number(tFinalHarmonicMotion));
    needAnalysis = true;

    numStepHarmonic = tFinalHarmonicMotion/dtHarmonicMotion+1;
    harmonicData = new Vector(numStepHarmonic);
    for (int i=0; i<numStepHarmonic; i++)
        (*harmonicData)(i) = magHarmonicMotion*sin(2*3.14159*i*dtHarmonicMotion/periodHarmonicMotion);

    this->reset();

    theNodeResponse->setItem(numF);
    theForceDispResponse->setItem(1);
    theForceTimeResponse->setItem(1);

    floorMassFrame->setVisible(false);
    storyPropertiesFrame->setVisible(false);
    spreadSheetFrame->setVisible(true);
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

void
MainWindow::on_magHarmonicChanged()
{
    QString newText =  magHarmonic->text();
    double  newDouble = newText.toDouble();
    if (newDouble != magHarmonicMotion) {
        magHarmonicMotion = newDouble;

        for (int i=0; i<numStepHarmonic; i++)
            (*harmonicData)(i) = magHarmonicMotion*sin(2*3.14159*i*dtHarmonicMotion/periodHarmonicMotion);

        this->setData(numStepHarmonic, dtHarmonicMotion, harmonicData);
        needAnalysis = true;
        this->reset();
 }
}

void
MainWindow::on_periodHarmonicChanged()
{
    QString newText =  periodHarmonic->text();
    double  newDouble = newText.toDouble();
    if (newDouble != periodHarmonicMotion) {
        periodHarmonicMotion = newDouble;
        needAnalysis = true;

        for (int i=0; i<numStepHarmonic; i++)
            (*harmonicData)(i) = magHarmonicMotion*sin(2*3.14159*i*dtHarmonicMotion/periodHarmonicMotion);

        this->setData(numStepHarmonic, dtHarmonicMotion, harmonicData);
    }
    this->reset();
}


void
MainWindow::on_dtHarmonicChanged()
{
    QString newText =  dtHarmonic->text();
    double  newDouble = newText.toDouble();
    if (newDouble != dtHarmonicMotion) {
        dtHarmonicMotion = newDouble;

        needAnalysis = true;
        if (harmonicData != 0)
            delete harmonicData;
        numStepHarmonic = tFinalHarmonicMotion/dtHarmonicMotion+1;
        harmonicData = new Vector(numStepHarmonic);
        for (int i=0; i<numStepHarmonic; i++)
            (*harmonicData)(i) = magHarmonicMotion*sin(2*3.14159*i*dtHarmonicMotion/periodHarmonicMotion);

        this->setData(numStepHarmonic, dtHarmonicMotion, harmonicData);
        this->reset();
 }
}


void
MainWindow::on_tFinalHarmonicChanged()
{
    QString newText =  tFinalHarmonic->text();
    double  newDouble = newText.toDouble();
    if (newDouble != tFinalHarmonicMotion) {
        tFinalHarmonicMotion = newDouble;
        if (harmonicData != 0)
            delete harmonicData;
        numStepHarmonic = tFinalHarmonicMotion/dtHarmonicMotion+1;
        harmonicData = new Vector(numStepHarmonic);
        for (int i=0; i<numStepHarmonic; i++)
            (*harmonicData)(i) = magHarmonicMotion*sin(2*3.14159*i*dtHarmonicMotion/periodHarmonicMotion);

        this->setData(numStepHarmonic, dtHarmonicMotion, harmonicData);

        needAnalysis = true;
        this->reset();
 }
}


void MainWindow::on_inFloors_editingFinished()
{
    QString textFloors =  inFloors->text();
    int numFloorsText = textFloors.toInt();
    if (numFloorsText != numFloors) {
        this->setBasicModel(numFloorsText, buildingW, buildingH, storyK, dampingRatio, g);
        floorSelected = -1;
        storySelected = -1;
        fMinSelected = -1;
        fMaxSelected = -1;
        sMinSelected = -1;
        sMaxSelected = -1;
       // this->setSelectionBoundary(-1.,-1.);
    }
}

void MainWindow::on_inWeight_editingFinished()
{
    QString textW =  inWeight->text();
    double textToDoubleW = textW.toDouble();
    if (textToDoubleW != buildingW) {
        buildingW = textToDoubleW;
        double floorW = buildingW/(numFloors);

        for (int i=0; i<numFloors; i++) {
            weights[i] = floorW;
        }
         this->updatePeriod();
        needAnalysis = true;
        this->reset();

        //inHeight->setFocus();
    }
}

void MainWindow::on_inHeight_editingFinished()
{
    QString textH =  inHeight->text();
    if (textH.isNull())
        return;
    double textToDoubleH = textH.toDouble();
    if (textToDoubleH != buildingH) {
        buildingH = textToDoubleH;

        for (int i=0; i<numFloors; i++) {
            floorHeights[i] = i*buildingH/(1.*numFloors);
            storyHeights[i] = buildingH/(1.*numFloors);
        }
        floorHeights[numFloors] = buildingH;

        this->updatePeriod();
        needAnalysis = true;
        this->reset();
    }
}


void MainWindow::on_inK_editingFinished()
{
    QString text =  inK->text();
    if (text.isNull())
        return;

    double textToDouble = text.toDouble();
    if (textToDouble != storyK) {
        storyK = textToDouble;

        for (int i=0; i<numFloors; i++) {
            k[i] = storyK;
        }

        this->updatePeriod();
        needAnalysis = true;
        this->reset();
        //this->setBasicModel(numFloors, buildingW, buildingH, storyK, dampingRatio, g);
    }
}

void MainWindow::on_inDamping_editingFinished()
{
    QString text =  inDamping->text();
    if (text.isNull())
        return;

    double textToDouble = text.toDouble();
    if (dampingRatio != textToDouble) {
        if (textToDouble < 0 || textToDouble > 1.0) {
            inDamping->setText(QString::number(dampingRatio));
            return;
        }
        dampingRatio = textToDouble;
        //    inGravity->setFocus();
        for (int i=0; i<numFloors; i++) {
            dampRatios[i]=dampingRatio;
        }
        this->reset();
    }
}

void MainWindow::on_scaleFactor_editingFinished()
{
    QString text =  scaleFactorEQ->text();
    if (text.isNull())
        return;

    double textToDouble = text.toDouble();
    if (scaleFactor != textToDouble) {
        scaleFactor = textToDouble;
        theCurrentRecord->setScaleFactor(scaleFactor);
        needAnalysis=true;

        double maxValue=0;
        for (int i = 0; i < numSteps; ++i) {
            double value = (*motionData)[i] * scaleFactor;
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


        QString textText("pga: "); textText.append(QString::number(maxValue,'g',2)); textText.append(tr("g"));
        earthquakeText->setText(textText);

        earthquakePlot->replot();
        earthquakePlot->update();

        this->reset();
    }
}

void MainWindow::on_inGravity_editingFinished()
{

}

void MainWindow::on_inFloorWeight_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inFloorWeight->text();
    if (text.isNull())
        return;

    double textToDouble = text.toDouble();
    for (int i=fMinSelected; i<=fMaxSelected; i++)
        weights[i] = textToDouble;

    buildingW = 0;
    for (int i=0; i<numFloors; i++)
        buildingW = buildingW+weights[i];

    inWeight->setText(QString::number(buildingW));
    this->reset();
}



void MainWindow::on_inStoryHeight_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryHeight->text();
    if (text.isNull())
        return;

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



    bool needReset = false;
    for (int i=0; i<=numFloors; i++) {
        if (floorHeights[i] != newFloorHeights[i]){
            needReset = true;
            i=numFloors+1;
        }
    }


    delete [] floorHeights;
    floorHeights = newFloorHeights;

    // move focus, update graphic and set analysis flag
    //  inStoryK->setFocus();
    if (needReset == true) {
        // delete old array and reset pointer
        this->reset();
        // delete old array and reset pointer
        buildingH = newFloorHeights[numFloors];
        inHeight->setText(QString::number(buildingH));
    }
}

void MainWindow::on_inStoryK_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryK->text();
    if (text.isNull())
        return;

    double textToDouble = text.toDouble();
    bool needReset = false;
    for (int i=sMinSelected; i<=sMaxSelected; i++) {
        if (k[i] != textToDouble) {
            k[i] = textToDouble;
            needReset = true;
        }
    }

    //inStoryFy->setFocus();
    if (needReset == true)
        this->reset();
}

void MainWindow::on_inStoryFy_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryFy->text();
    if (text.isNull())
        return;

    bool needReset = false;
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++) {
        if (fabs(fy[i]-textToDouble) > 1e-12) {
            fy[i] = textToDouble;
            needReset = true;
        }
    }
    //inStoryB->setFocus();
    if (needReset == true)
        this->reset();
}

void MainWindow::on_inStoryB_editingFinished()
{
    if (updatingPropertiesTable == true)
        return;

    QString text =  inStoryB->text();
    if (text.isNull())
        return;

    bool needReset = false;
    double textToDouble = text.toDouble();
    for (int i=sMinSelected; i<=sMaxSelected; i++)
        if (b[i] != textToDouble) {
            b[i] = textToDouble;
            needReset = true;
        }

    //inStoryHeight->setFocus();
    if (needReset == true)
        this->reset();
}


void MainWindow::doAnalysis()
{ 
    if (needAnalysis == true && analysisFailed == false) {

      //        qDebug() << "doANALYSIS";

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
        ZeroLength **theElements = new ZeroLength *[numFloors];
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
            theElements[i-1] = theEle;
            theDomain.addElement(theEle);
            delete theMat; // each ele makes it's own copy
        }

        //
        // create load pattern and add loads
        //
        PathSeries *theSeries;
        theSeries = new PathSeries(1, *motionData, dt, g*scaleFactor);
          /*
        if (motionTypeValue == 0) {
            theSeries = new PathSeries(1, *eqData, dt, g*scaleFactor);
        } else {
            theSeries = new PathSeries(1, *harmonicData, dtHarmonicMotion, g*magHarmonicMotion);
        }
        */
        GroundMotion *theGroundMotion = new GroundMotion(0,0,theSeries);
        LoadPattern *theLoadPattern = new UniformExcitation(*theGroundMotion, 0, 1);

        //   theLoadPattern->setTimeSeries(theTimeSeries);
        //   static Vector load(1); load.Zero(); load(0) = 1;
        //   NodalLoad *theLoad = new NodalLoad(0, numFloors, load);
        //   theLoadPattern->addNodalLoad(theLoad);

        theDomain.addLoadPattern(theLoadPattern);

        //
        // create the analysis
        //

        AnalysisModel     *theModel = new AnalysisModel();
        CTestNormDispIncr *theTest = new CTestNormDispIncr(1.0e-10, 20, 0);
        EquiSolnAlgo      *theSolnAlgo = new NewtonRaphson();//INITIAL_TANGENT);
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
        if (eigValues->Size() != numFloors) {
                eigValues->resize(numFloors);
        }
         *eigValues = theEig;
        if (ok == 0)
            for (int i=0; i<numFloors; i++)
                if (theEig(i) <= 0)
                    ok = -1;

        int numCombo = periodComboBox->count();

        if (numCombo != numFloors) {
            periodComboBox->clear();
            QString t1 = QString("Fundamental Period");
            periodComboBox->addItem(t1);
            for (int i=1; i<numFloors; i++) {
                QString t1 = QString("T") + QString::number(i+1);
                periodComboBox->addItem(t1);
            }
        }

        if (ok != 0) {
            needAnalysis = false;
            analysisFailed = true;
            for (int i=0; i<numSteps; i++) {
                for (int j=0; j<numFloors+1; j++) {
                    dispResponses[j][i] = 0;
                    if (j < numFloors) {
                        storyForceResponses[j][i]=0;
                        storyDriftResponses[j][i]=0;
                    }
                }
            }
            maxDispLabel->setText(QString().setNum(0.0,'f',2));
            currentPeriod->setText(QString(tr("undefined")));
            delete [] theNodes;
            delete [] theElements;

            QMessageBox::warning(this, tr("Application"),
                                 tr("Eigenvalue Analysis Failed.<p> Possible Causes: negstive mass, negative story stiffness, or "
                                    "if PDelta is included, a story stiffness - axial load divided by L is negtive resulting "
                                    "in non postive definite stiffness matrix"));

            return;
        }

        Vector dampValues(numFloors);
        for (int i=0; i<numFloors; i++) {
            dampValues(i)=dampRatios[i];
        }
        theDomain.setModalDampingFactors(&dampValues);
        double T1 = 2*3.14159/sqrt(theEig(0));

        maxDisp = 0;
        for (int i=0; i<=numSteps; i++) { // <= due to adding 0 at start
            int ok = theAnalysis.analyze(1, dt);
            if (ok != 0) {
                needAnalysis = false;
                analysisFailed = true;
                for (int k=i; k<numSteps; k++) {
                    for (int j=0; j<numFloors+1; j++) {
                        dispResponses[j][i] = 0;
                        if (j < numFloors) {
                            storyForceResponses[j][k]=0;
                            storyDriftResponses[j][k]=0;
                        }
                    }

                }
                QMessageBox::warning(this, tr("Application"),
                                     tr("Transient Analysis Failed"));
                break;
            }
            for (int j=0; j<numFloors+1; j++) {
                double nodeDisp = theNodes[j]->getDisp()(0);
                dispResponses[j][i] = nodeDisp;
                if (fabs(nodeDisp) > maxDisp)
                    maxDisp = fabs(nodeDisp);

                if (j < numFloors) {

                    storyForceResponses[j][i]= theElements[j]->getForce();
                    storyDriftResponses[j][i] = theElements[j]->getDrift();
                }
            }
        }
        // clean up memory
        delete [] theNodes;
        delete [] theElements;

        // reset values, i.e. slider position, current displayed step...
        maxDispLabel->setText(QString().setNum(maxDisp,'f',2));

        int num = periodComboBox->currentIndex();
        double eigenValue = (*eigValues)(num);
        if (eigenValue <= 0) {
            currentPeriod->setText(QString("undefined"));
        } else {
            double period = 2*3.14159/sqrt((eigenValue));
            currentPeriod->setText(QString().setNum(period,'f',2));
        }
       // currentPeriod->setText(QString().setNum(T1,'f',2));
        needAnalysis = false;
        currentStep = 0;
        //  groupTracer->setGraphKey(0);
        //      slider->setSliderPosition(0);
        myGL->update();

        analysisFailed = false;
        int nodeResponseFloor = theNodeResponse->getItem();
        int storyForceTime = theForceTimeResponse->getItem() -1;
        //int storyForceDrift = theForceDispResponse->getItem() -1;

        nodeResponseValues.resize(numSteps);
        storyForceValues.resize(numSteps);
        storyDriftValues.resize(numSteps);

        for (int i = 0; i < numSteps; ++i) {
            nodeResponseValues[i]=dispResponses[nodeResponseFloor][i];
            storyForceValues[i]=storyForceResponses[storyForceTime][i];
            storyDriftValues[i]=storyDriftResponses[storyForceTime][i];
            // qDebug() << i*dt << " " << storyForceResponses[storyForceTime][i] << " " << storyDriftResponses[storyForceTime][i];
        }


        theNodeResponse->setData(nodeResponseValues,time,numSteps,dt);
        theForceTimeResponse->setData(storyForceValues,time,numSteps,dt);
        theForceDispResponse->setData(storyForceValues,storyDriftValues,numSteps);
    }
}

void
MainWindow::setResponse(int floor, int mainItem)
{
    if (mainItem == 0) {
        if (floor > 0 && floor <= numFloors) {
            for (int i = 0; i < numSteps; ++i) {
                nodeResponseValues[i]=dispResponses[floor][i];
            }
            theNodeResponse->setData(nodeResponseValues,time,numSteps,dt);
        }
    } else if (mainItem == 1 || mainItem == 2) {
        if (floor > 0 && floor <= numFloors) {
            for (int i = 0; i < numSteps; ++i) {
                storyForceValues[i]=storyForceResponses[floor-1][i];
                storyDriftValues[i]=storyDriftResponses[floor-1][i];
            }
            theForceTimeResponse->setData(storyForceValues,time,numSteps,dt);
            theForceDispResponse->setData(storyForceValues,storyDriftValues,numSteps);
            if (mainItem == 1)
                theForceDispResponse->setItem(floor);
            else
                theForceTimeResponse->setItem(floor);
        }
    }
}

void MainWindow::reset() {

    analysisFailed = false;
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
            if (weights[row] == textToDouble)
                return;

            weights[row] = textToDouble;
            buildingW = 0;
            for (int i=0; i<numFloors; i++)
                buildingW += weights[i];
            inWeight->setText(QString::number(buildingW));
        } else if (column  == 1) {
            if (storyHeights[row] == textToDouble)
                return;

            storyHeights[row] = textToDouble;
            for (int i=row; i<numFloors; i++)
                floorHeights[i+1] = floorHeights[i]+storyHeights[i];
            buildingH = floorHeights[numFloors];
            inHeight->setText(QString::number(buildingH));

        } else if (column == 2) {
            if (k[row] == textToDouble)
                return;

            k[row] = textToDouble;
        } else if (column == 3) {         
            if (fy[row] == textToDouble)
                return;

            fy[row] = textToDouble;
        } else if (column == 4) {
            if (b[row] == textToDouble)
                return;

            b[row] = textToDouble;
        } else {
            if (dampRatios[row] == textToDouble)
                return;

            dampRatios[row] = textToDouble;
        }


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
        fMinSelected = -1;
        fMaxSelected = -1;
        sMinSelected = -1;
        sMaxSelected = -1;


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
        spreadSheetFrame->setVisible(true);
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
    return;
}

void MainWindow::on_PeriodSelectionChanged(const QString &arg1) {
    int num = periodComboBox->currentIndex();
    double eigenValue = (*eigValues)(num);
    if (eigenValue <= 0) {
     currentPeriod->setText(QString("undefined"));
    } else {
    double period = 2*3.14159/sqrt((eigenValue));
    currentPeriod->setText(QString().setNum(period,'f',2));
    }
}


void MainWindow::on_motionTypeSelectionChanged(const QString &arg1)
{
    if (arg1 == QString(tr("Harmonic Motion"))) {
        motionTypeValue = 1;
        eqMotionFrame->setVisible(false);
        harmonicMotionFrame->setVisible(true);
        scaleFactor=1.0;
        this->setData(numStepHarmonic, dtHarmonicMotion, harmonicData);

    } else {
        motionTypeValue = 0;
        harmonicMotionFrame->setVisible(false);
        eqMotionFrame->setVisible(true);
        scaleFactor=(scaleFactorEQ->text()).toDouble();

         this->setData(numStepEarthquake, dtEarthquakeMotion, eqData);
    }
   // this->doAnalysis();
    needAnalysis = true;
    this->reset();
    return;
}


void MainWindow::setData(int nStep, double deltaT, Vector *data) {

    numSteps = nStep;
    dt = deltaT;
    motionData = data;

    double maxValue = 0;

    if (dispResponses != 0) {
        for (int j=0; j<numFloors+1; j++)
            delete [] dispResponses[j];
        delete [] dispResponses;
    }
    if (storyForceResponses != 0) {
        for (int j=0; j<numFloors; j++)
            delete [] storyForceResponses[j];
        delete [] storyForceResponses;
    }
    if (storyDriftResponses != 0) {
        for (int j=0; j<numFloors; j++)
            delete [] storyDriftResponses[j];
        delete [] storyDriftResponses;
    }

    dispResponses = new double *[numFloors+1];
    storyForceResponses = new double *[numFloors];
    storyDriftResponses = new double *[numFloors];

    for (int i=0; i<numFloors+1; i++) {
        dispResponses[i] = new double[numSteps+1]; // +1 as doing 0 at start
        if (i<numFloors) {
            storyForceResponses[i] = new double[numSteps+1];
            storyDriftResponses[i] = new double[numSteps+1];
        }
    }


    excitationValues.resize(numSteps);
    time.resize(numSteps);

    for (int i = 0; i < numSteps; ++i) {
        double value = (*motionData)[i] * scaleFactor;
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


    QString textText("pga: "); textText.append(QString::number(maxValue,'g',2)); textText.append(tr("g"));
    earthquakeText->setText(textText);

    earthquakePlot->replot();
    earthquakePlot->update();


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

    needAnalysis = true;
  //  this->reset();
  //  myGL->update();
}


void MainWindow::on_inEarthquakeMotionSelectionChanged(const QString &arg1)
{
    std::map<QString, EarthquakeRecord *>::iterator it;
    it = records.find(arg1);
    if (it != records.end()) {
        theCurrentRecord = records.at(arg1);
        numStepEarthquake =  theCurrentRecord->numSteps;
        dtEarthquakeMotion =  theCurrentRecord->dt;
        eqData = theCurrentRecord->data;
        scaleFactor=theCurrentRecord->getScaleFactor();
        scaleFactorEQ->setText(QString::number(scaleFactor));
        this->setData(numStepEarthquake, dtEarthquakeMotion, eqData);
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
    this->setCurrentFile(fileName);
}

void MainWindow::resetFile()
{
    // reset to original
    this->setBasicModel(5, 5*100, 5*144, 31.54, .05, 386.4);
    this->setSelectionBoundary(-1.,-1.);

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
    if (theValue.isNull()) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot find \"name\" attribute in file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(inputMotionName),
                                  file.errorString()));
        return;

    }
    QString name =theValue.toString();

    // create a record and add to records map
    EarthquakeRecord *theRecord = new EarthquakeRecord();
    int ok = theRecord->inputFromJSON(jsonObject);

    if (ok == 0) {
        // inser into records list
        records.insert(std::make_pair(name, theRecord));

        // add the motion to ComboBox, & set it current
        eqMotion->addItem(name);
        int index = eqMotion->findText(name);
        eqMotion->setCurrentIndex(index);

        needAnalysis = true;
        analysisFailed = false;

    } else {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot find an attribute needed in file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(inputMotionName),
                                  file.errorString()));
    }

    // close file & return
    file.close();
    return;
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
    json["currentMotion"]=eqMotion->currentText();
    json["currentMotionIndex"]=eqMotion->currentIndex();

    json["periodHarmonic"]=periodHarmonicMotion;
    json["magHarmonic"]=magHarmonicMotion;
    json["tFinalHarmonic"]=tFinalHarmonicMotion;
    json["dtHarmonic"]=dtHarmonicMotion;

    json["motionType"]=motionType->currentIndex();

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

    QJsonArray motionsArray;
    int numMotions = eqMotion->count();
    std::map<QString, EarthquakeRecord *>::iterator iter;
    for (int i=0; i<numMotions; i++) {
        QString eqName = eqMotion->itemText(i);
        iter = records.find(eqName);
        if (iter != records.end()) {
            QJsonObject obj;
            EarthquakeRecord *theRecord = iter->second ;
            theRecord->outputToJSON(obj);
            motionsArray.append(obj);
        }
    }

    json["records"]=motionsArray;

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
        <p>\
        The source code is licensed under a BSD 2-Clause License:<p>\
        \"Copyright (c) 2017-2018, The Regents of the University of California (Regents).\"\
        All rights reserved.<p>\
        <p>\
        Redistribution and use in source and binary forms, with or without \
        modification, are permitted provided that the following conditions are met:\
        <p>\
         1. Redistributions of source code must retain the above copyright notice, this\
         list of conditions and the following disclaimer.\
         \
         \
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
            The compiled binary form of this application is licensed under a GPL Version 3 license.\
            The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
            included in the packaging of this application. \
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            This software makes use of the QT packages (unmodified): core, gui, widgets and network\
                                                                     <p>\
                                                                     QT is copyright \"The Qt Company Ltd&quot; and licensed under the GNU Lesser General \
                                                                     Public License (version 3) which references the GNU General Public License (version 3)\
      <p>\
      The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
      included in the packaging of this application. \
      <p>\
      ------------------------------------------------------------------------------------\
      <p>\
      This software makes use of the OpenSees Software Framework. OpenSees is copyright \"The Regents of the University of \
      California\". OpenSees is open-source software whose license can be\
      found at http://opensees.berkeley.edu.\
      <p>\
      ";


         QMessageBox msgBox;
    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(textCopyright);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
    msgBox.exec();

}


void MainWindow::version()
{
    QMessageBox::about(this, tr("Version"),
                       tr("Version 1.1 "));
}

void MainWindow::about()
{
    QString textAbout = "\
            This is the Multiple Degree of Freedom (MDOF) tool.  It allows the user to explore the effects of\
            different building parameters and ground motions on the time history response of a building. <p> \
            The building is represented by a shear building model: an idealization of a structure in which the mass \
            is lumped at the floor levels and the beams are assumed infinitely stiff in flexure and axially inextensible,\
            and the columns are axially inextensible.  The user inputs the floor weights and story properties (stiffness, \
                                                                                                               yield strength, hardening ratio) of the stories, and a damping ratio for the structure. Individual floor and \
            story values are possible by user selecting an an appropriate area in the graphic around area of interest.\
            In addition nonlinear effects due to P-Delta and soft story mechanisms can be studied.\
            <p>\
            All units are in sec, kips, inches.\
            <p>\
            For this application the equations of motions are set up using the uniform excitation approach, \
            i.e. MA + CV + KU = -MAg. These equations are solved using the Newmark constant acceleration method and \
            Newton-Raphson solution algorithm.  <p>\
            Additional motions can be added by user. The units for these additional motions must be in g. An\
            example is provided at https://github.com/NHERI-SimCenter/MDOF/blob/master/example/elCentro.json\
            <p>\
            This tool does not stop you the user from inputting values that will cause the analysis to fail. \
            If the analysis fails, a warning message will appear and the tool will not revert back to a working \
            set of parameters.\
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
   // QDesktopServices::openUrl(QUrl("https://github.com/NHERI-SimCenter/MDOF/issues", QUrl::TolerantMode));
 QDesktopServices::openUrl(QUrl("https://www.designsafe-ci.org/help/new-ticket/", QUrl::TolerantMode));
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

  //
  // clean up old
  //

  if (dispResponses != 0) {
    for (int j=0; j<numFloors+1; j++)
      delete [] dispResponses[j];
    delete [] dispResponses;
  }
  if (storyForceResponses != 0) {
    for (int j=0; j<numFloors; j++)
      delete [] storyForceResponses[j];
    delete [] storyForceResponses;
  }
  if (storyDriftResponses != 0) {
    for (int j=0; j<numFloors; j++)
      delete [] storyDriftResponses[j];
    delete [] storyDriftResponses;
        }
  
  dispResponses = 0;
  storyForceResponses = 0;
  storyDriftResponses = 0;

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
    //inGravity->setText(QString::number(g));

    theValue = jsonObject["dampingRatio"];
    dampingRatio=theValue.toDouble();
    inDamping->setText(QString::number(dampingRatio));

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


    //
    // clear records and inMotion combo box
    //

    std::map<QString, EarthquakeRecord *>::iterator iter;

    // delete the earthqkaes before clear as clear does not invoke destructor if pointers
    for (iter = records.begin(); iter != records.end(); ++iter ) {
        delete iter->second;
    }

    records.clear();
    eqMotion->clear();

    theValue = jsonObject["records"];
    theArray=theValue.toArray();

    //
    // now read records from file and populate records and comoboBox with new
    //

    for (int i=0; i<theArray.size(); i++) {
        EarthquakeRecord *theRecord = new EarthquakeRecord();
        QJsonObject theEarthquakeObj = theArray.at(i).toObject();
        theRecord->inputFromJSON(theEarthquakeObj);
        records.insert(std::make_pair(theRecord->name, theRecord));
        eqMotion->addItem(theRecord->name);
    }

    //
    // set current record to what was saved in file
    //

    theValue = jsonObject["currentMotionIndex"];
    eqMotion->setCurrentIndex(theValue.toInt());
    theValue = jsonObject["currentMotion"];


    //
    // read harmonic data
    //

    theValue = jsonObject["periodHarmonic"];
    periodHarmonicMotion=theValue.toDouble();
    periodHarmonic->setText(QString::number(periodHarmonicMotion));

    theValue = jsonObject["magHarmonic"];
    magHarmonicMotion=theValue.toDouble();
    magHarmonic->setText(QString::number(magHarmonicMotion));

    theValue = jsonObject["tFinalHarmonic"];
    tFinalHarmonicMotion=theValue.toDouble();
    tFinalHarmonic->setText(QString::number(tFinalHarmonicMotion));

    theValue = jsonObject["dtHarmonic"];
    dtHarmonicMotion=theValue.toDouble();
    dtHarmonic->setText(QString::number(dtHarmonicMotion));


   // json["motionType"]=motionType->currentIndex();
    theValue = jsonObject["motionType"];
    motionType->setCurrentIndex(theValue.toInt());
    //dtHarmonicMotion=theValue.toDouble();
    //dtHarmonic->setText(QString::number(dtHarmonicMotion));


    this->reset();
    this->on_inEarthquakeMotionSelectionChanged(theValue.toString());
    theNodeResponse->setItem(numFloors);
    theForceDispResponse->setItem(1);
    theForceTimeResponse->setItem(1);

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

    QAction *newAction = new QAction(tr("&Reset"), this);
    newAction->setShortcuts(QKeySequence::New);
    newAction->setStatusTip(tr("Create a new file"));
    connect(newAction, &QAction::triggered, this, &MainWindow::resetFile);
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
    connect(saveAsAction, &QAction::triggered, this, &MainWindow::saveAs);
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

    QString tLabel("Time");
    QString dLabel("Relative Displacement");
    QString fLabel("Floor");
    QString sLabel("Story");
    QString forceLabel("Shear Force");
    QString dispLabel("Displacement");

    theNodeResponse = new ResponseWidget(this, 0, fLabel, tLabel, dLabel);
    QDockWidget *nodeResponseDock = new QDockWidget(tr("Floor Displacement History"), this);
    nodeResponseDock->setWidget(theNodeResponse);
    nodeResponseDock->setAllowedAreas(Qt::NoDockWidgetArea);
    nodeResponseDock->setFloating(true);
    nodeResponseDock->close();
    viewMenu->addAction(nodeResponseDock->toggleViewAction());

    theForceTimeResponse = new ResponseWidget(this, 1, sLabel, tLabel, forceLabel);
    QDockWidget *forceTimeResponseDock = new QDockWidget(tr("Story Force History"), this);
    forceTimeResponseDock->setWidget(theForceTimeResponse);
    forceTimeResponseDock->setAllowedAreas(Qt::NoDockWidgetArea);
    forceTimeResponseDock->setFloating(true);
    forceTimeResponseDock->close();
    viewMenu->addAction(forceTimeResponseDock->toggleViewAction());


    theForceDispResponse = new ResponseWidget(this, 2, sLabel, dispLabel,forceLabel);
    QDockWidget *forceDriftResponseDock = new QDockWidget(tr("Story Force-Displacement"), this);
    forceDriftResponseDock->setWidget(theForceDispResponse);
    forceDriftResponseDock->setAllowedAreas(Qt::NoDockWidgetArea);
    forceDriftResponseDock->setFloating(true);
    forceDriftResponseDock->close();
    viewMenu->addAction(forceDriftResponseDock->toggleViewAction());


    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
    QAction *infoAct = helpMenu->addAction(tr("&About"), this, &MainWindow::about);
    QAction *submitAct = helpMenu->addAction(tr("&Provide Feedback"), this, &MainWindow::submitFeedback);
    //aboutAct->setStatusTip(tr("Show the application's About box"));
    QAction *aboutAct = helpMenu->addAction(tr("&Version"), this, &MainWindow::version);
    //aboutAct->setStatusTip(tr("Show the application's About box"));


    QAction *copyrightAct = helpMenu->addAction(tr("&License"), this, &MainWindow::copyright);
    //aboutAct->setStatusTip(tr("Show the application's About box"));

}

void MainWindow::viewNodeResponse(){

}

void MainWindow::viewStoryResponse(){

}

void MainWindow::createHeaderBox() {

    HeaderWidget *header = new HeaderWidget();
    header->setHeadingText(tr("Multiple Degrees of Freedom Application"));

    largeLayout->addWidget(header);
}

void MainWindow::createFooterBox() {

    FooterWidget *footer = new FooterWidget();

    largeLayout->addWidget(footer);
}

void MainWindow::createInputPanel() {
    inputLayout = new QVBoxLayout;

    //
    // some QStrings to avoid duplication of units
    //

    QString blank(tr("   "));
    QString kips(tr("k"  ));
    QString g(tr("g"  ));
    QString kipsInch(tr("k/in"));
    QString inch(tr("in  "));
    QString sec(tr("sec"));
    QString percent(tr("\%   "));

    //
    // Create a section line + title + add
    // styleSheet

    QHBoxLayout *inputMotionType = new QHBoxLayout();
    SectionTitle *inTitle = new SectionTitle(this);
    inTitle->setTitle(tr("Input Motion"));
    motionType = new QComboBox();
    inputMotionType->addWidget(inTitle);
    inputMotionType->addStretch();
    inputMotionType->addWidget(motionType);
    QString eqString("Earthquake Motion");
    QString harmonicString("Harmonic Motion");
    motionType->addItem(eqString);
    motionType->addItem(harmonicString);

    inputLayout->addLayout(inputMotionType);
    motionTypeValue = 0;

    // JUST EARTHQUAKE inputLayout->addWidget(inTitle);

    //
    // create the frame for the earthquake motion selection
    //

    eqMotionFrame = new QFrame(); //styleSheet
    eqMotionFrame->setObjectName(QString::fromUtf8("inputMotion")); //styleSheet

    QGridLayout *inputMotionLayout = new QGridLayout();
    //QHBoxLayout *inputMotionLayout = new QHBoxLayout();
    QLabel *sectionTitle = new QLabel();
    sectionTitle->setText(tr("Input Motion Title"));
    sectionTitle->setObjectName(QString::fromUtf8("sectionTitle")); //styleSheet
    QLabel *entryLabel = new QLabel();
    entryLabel->setText(tr("Input Motion"));

    eqMotion = new QComboBox();

    QLabel *scaleLabel = new QLabel(tr("Scale Factor"));
    scaleFactorEQ = new QLineEdit();
    scaleFactorEQ->setMaximumWidth(100);

    inputMotionLayout->addWidget(entryLabel,0,0);
    inputMotionLayout->addWidget(eqMotion,0,1);
    inputMotionLayout->addWidget(scaleLabel,0,3);
    inputMotionLayout->addWidget(scaleFactorEQ,0,4);
   // inputMotionLayout->addStretch();


    addMotion = new QPushButton("Add");
    inputMotionLayout->addWidget(addMotion,1,4);
    eqMotionFrame->setLayout(inputMotionLayout);
    inputMotionLayout->setColumnStretch(2,1);
    // inputMotion->setFrameStyle(QFrame::Raised);
    eqMotionFrame->setLineWidth(1);
    eqMotionFrame->setFrameShape(QFrame::Box);
    eqMotionFrame->setVisible(true);

    inputLayout->addWidget(eqMotionFrame);

    //
    // create the frame for the harmonic motion
    //

    // frame to hold it with QGridLayout
    harmonicMotionFrame = new QFrame(); //style sheet
    harmonicMotionFrame->setObjectName(QString::fromUtf8("inputMotion")); //styleSheet

    QGridLayout *harmonicLayout = new QGridLayout();

    // the widgets for the layout
    QLabel *periodLabel = new QLabel(tr("Period:"));
    QLabel *periodUnit = new QLabel(sec);
    periodHarmonic = new QLineEdit();

    QLabel *magLabel = new QLabel(tr("PGA:"));
    QLabel *magUnit = new QLabel(g);
    magHarmonic = new QLineEdit();

    QLabel *dTLabel = new QLabel(tr("delta T"));
    dtHarmonic = new QLineEdit();
    QLabel *dTUnit = new QLabel(sec);

    QLabel *tFinalLabel = new QLabel(tr("tFinal"));
    tFinalHarmonic = new QLineEdit();
    QLabel *tFinalUnit = new QLabel(sec);

    // add the widgets to the layout
    harmonicLayout->addWidget(periodLabel,0,0);
    harmonicLayout->addWidget(periodHarmonic,0,1);
    harmonicLayout->addWidget(periodUnit,0,2);

    harmonicLayout->addWidget(magLabel,1,0);
    harmonicLayout->addWidget(magHarmonic,1,1);
    harmonicLayout->addWidget(magUnit,1,2);

    harmonicLayout->setColumnStretch(3,1);

    harmonicLayout->addWidget(dTLabel,0,4);
    harmonicLayout->addWidget(dtHarmonic,0,5);
    harmonicLayout->addWidget(dTUnit,0,6);

    harmonicLayout->addWidget(tFinalLabel,1,4);
    harmonicLayout->addWidget(tFinalHarmonic,1,5);
    harmonicLayout->addWidget(tFinalUnit,1,6);

    harmonicMotionFrame->setVisible(false);
    harmonicMotionFrame->setLayout(harmonicLayout);
    harmonicMotionFrame->setLineWidth(1);
    harmonicMotionFrame->setFrameShape(QFrame::Box);

    inputLayout->addWidget(harmonicMotionFrame);

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

    SectionTitle *inBuilding = new SectionTitle(this);
    inBuilding->setTitle(tr("Buiding Properties"));
    inputLayout->addWidget(inBuilding);

    //
    // create to hold major model inputs
    //

    QFrame *mainProperties = new QFrame();
    mainProperties->setObjectName(QString::fromUtf8("mainProperties")); //styleSheet
    QVBoxLayout *mainPropertiesLayout = new QVBoxLayout();
    inFloors = createTextEntry(tr("Number Floors"), mainPropertiesLayout, 100, 100, &blank);
    inWeight = createTextEntry(tr("Building Weight"), mainPropertiesLayout, 100, 100, &kips);
    inHeight = createTextEntry(tr("Building Height"), mainPropertiesLayout, 100, 100, &inch);
    inK = createTextEntry(tr("Story Stiffness"), mainPropertiesLayout, 100, 100, &kipsInch);
    inDamping = createTextEntry(tr("Damping Ratio"), mainPropertiesLayout, 100, 100, &percent);
    //inGravity =  createTextEntry(tr("Gravity"), mainPropertiesLayout);
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
    spreadSheetFrame->setObjectName(QString::fromUtf8("inputMotion")); //styleSheet
    spreadSheetFrame->setLayout(spreadsheetFrameLayout);
    spreadSheetFrame->setLineWidth(1);
    spreadSheetFrame->setFrameShape(QFrame::Box);

    inputLayout->addWidget(spreadSheetFrame,1);
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

    scaleFactorEQ->setValidator(new QDoubleValidator);

    periodHarmonic->setValidator(new QDoubleValidator);
    magHarmonic->setValidator(new QDoubleValidator);
    dtHarmonic->setValidator(new QDoubleValidator);
    tFinalHarmonic->setValidator(new QDoubleValidator);

    inFloors->setValidator(new QIntValidator);
    inWeight->setValidator(new QDoubleValidator);
    inHeight->setValidator(new QDoubleValidator);
    inK->setValidator(new QDoubleValidator);
    inDamping->setValidator(new QDoubleValidator);
    //inGravity->setValidator(new QDoubleValidator);
    inFloorWeight->setValidator(new QDoubleValidator);
    inStoryB->setValidator(new QDoubleValidator);
    inStoryFy->setValidator(new QDoubleValidator);
    inStoryHeight->setValidator(new QDoubleValidator);
    inStoryK->setValidator(new QDoubleValidator);

    //
    // connect signals & slots
    //

    connect(motionType, SIGNAL(currentIndexChanged(QString)), this, SLOT(on_motionTypeSelectionChanged(QString)));
    connect(scaleFactorEQ,SIGNAL(editingFinished()),this,SLOT(on_scaleFactor_editingFinished()));
    connect(magHarmonic,SIGNAL(editingFinished()), this, SLOT(on_magHarmonicChanged()));
    connect(periodHarmonic,SIGNAL(editingFinished()), this, SLOT(on_periodHarmonicChanged()));
    connect(dtHarmonic,SIGNAL(editingFinished()), this, SLOT(on_dtHarmonicChanged()));
    connect(tFinalHarmonic,SIGNAL(editingFinished()), this, SLOT(on_tFinalHarmonicChanged()));

    connect(pDeltaBox, SIGNAL(stateChanged(int)), this, SLOT(on_includePDeltaChanged(int)));
    connect(addMotion,SIGNAL(clicked()), this, SLOT(on_addMotion_clicked()));
    connect(inFloors,SIGNAL(editingFinished()), this, SLOT(on_inFloors_editingFinished()));
    connect(inWeight,SIGNAL(editingFinished()), this, SLOT(on_inWeight_editingFinished()));
    //connect(inWeight,SIGNAL(returnPressed()),this,SLOT(on_inWeight_editingFinished()));
    connect(inHeight,SIGNAL(editingFinished()), this, SLOT(on_inHeight_editingFinished()));
    connect(inK,SIGNAL(editingFinished()), this, SLOT(on_inK_editingFinished()));
    connect(inDamping,SIGNAL(editingFinished()), this, SLOT(on_inDamping_editingFinished()));
    connect(inFloorWeight,SIGNAL(editingFinished()), this, SLOT(on_inFloorWeight_editingFinished()));
    connect(inStoryB, SIGNAL(editingFinished()),this,SLOT(on_inStoryB_editingFinished()));
    connect(inStoryFy,SIGNAL(editingFinished()),this, SLOT(on_inStoryFy_editingFinished()));
    connect(inStoryHeight, SIGNAL(editingFinished()), this, SLOT(on_inStoryHeight_editingFinished()));
    connect(inStoryK, SIGNAL(editingFinished()), this, SLOT(on_inStoryK_editingFinished()));

    connect(eqMotion, SIGNAL(currentIndexChanged(QString)), this, SLOT(on_inEarthquakeMotionSelectionChanged(QString)));

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

    SectionTitle *outTitle = new SectionTitle(this);
    outTitle->setTitle(tr("Output"));
    outputLayout->addWidget(outTitle);


    QString inch(tr("in  "));
    QString sec(tr("sec"));


    // frame for basic outputs,
    //QFrame *outputMaxFrame = new QFrame();

    // frame for max disp / current period
    QFrame *firstOutput = new QFrame(); //styleSheet
    firstOutput->setObjectName(QString::fromUtf8("firstOutput"));
    QVBoxLayout *firstOutputLayout = new QVBoxLayout();
    maxDispLabel = createLabelEntry(tr("Max Disp"), firstOutputLayout, 100,100,&inch); //styleSheet
 //   currentPeriod= createLabelEntry(tr("Fundamental Period"),firstOutputLayout, 100,100,&sec); //styleSheet

    // create combobox for period selection
    QHBoxLayout *periodLayout = new QHBoxLayout();
    periodComboBox = new QComboBox();
    QString t1("FundamentalPeriod");
    periodComboBox->addItem(t1);

    QLabel *res = new QLabel();
    res->setMinimumWidth(100);
    res->setMaximumWidth(100);
    res->setAlignment(Qt::AlignRight);
    currentPeriod=res;

    periodLayout->addWidget(periodComboBox);
    periodLayout->addStretch();
    periodLayout->addWidget(currentPeriod);

    QLabel *unitLabel = new QLabel();
    unitLabel->setText(sec);
    unitLabel->setMinimumWidth(40);
    unitLabel->setMaximumWidth(100);
    periodLayout->addWidget(unitLabel);


    periodLayout->setSpacing(10);
    periodLayout->setMargin(0);

    firstOutputLayout->addLayout(periodLayout);

    firstOutput->setLayout(firstOutputLayout);
    firstOutput->setLineWidth(1);
    firstOutput->setFrameShape(QFrame::Box);
    outputLayout->addWidget(firstOutput);


    QVBoxLayout *outputMaxLayout = new QVBoxLayout();
    QLabel *vizTitle = new QLabel(); //styleSheet
    vizTitle->setText(tr("Visualization Section Title")); //styleSheet
    vizTitle->setObjectName(QString::fromUtf8("vizTitle")); //styleSheet

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
    myGL->setMinimumHeight(300);
    myGL->setMinimumWidth(250);
    myGL->setModel(this);
    outputLayout->addWidget(myGL,1.0);

    // input acceleration plot
    earthquakePlot=new QCustomPlot();
    earthquakePlot->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Minimum);
    earthquakePlot->setMinimumHeight(100);
    earthquakePlot->setMaximumHeight(100);

    earthquakeText = new QCPItemText(earthquakePlot);
    earthquakeText->position->setType(QCPItemPosition::ptAxisRectRatio);
    earthquakeText->position->setCoords(0.9,0.1);

    outputLayout->addWidget(earthquakePlot);

    // slider for manual movement
    slider=new QSlider(Qt::Horizontal);
    outputLayout->addWidget(slider);

    // output frame to show current time
    QFrame *outputDataFrame = new QFrame();
    outputDataFrame->setObjectName(QString::fromUtf8("outputDataFrame"));
    QVBoxLayout *outputDataLayout = new QVBoxLayout();
    currentTime = createLabelEntry(tr("Current Time"), outputDataLayout, 100,100, &sec);
    currentDisp = createLabelEntry(tr("Current Roof Disp"), outputDataLayout, 100,100,  &inch);
    outputDataFrame->setLayout(outputDataLayout);
    outputDataFrame->setLineWidth(1);
    outputDataFrame->setFrameShape(QFrame::Box);
    outputDataFrame->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    outputLayout->addWidget(outputDataFrame);

    // add layout to mainLayout and to largeLayout
    mainLayout->addLayout(outputLayout);
    largeLayout->addLayout(mainLayout); //styleSheet

    // signal and slot connects for slider
    connect(periodComboBox, SIGNAL(currentIndexChanged(QString)), this, SLOT(on_PeriodSelectionChanged(QString)));
    connect(slider, SIGNAL(sliderPressed()),  this, SLOT(on_slider_sliderPressed()));
    connect(slider, SIGNAL(sliderReleased()), this, SLOT(on_slider_sliderReleased()));
    connect(slider, SIGNAL(valueChanged(int)),this, SLOT(on_slider_valueChanged(int)));

    //outputLayout->addStretch();
}

