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

#include "MyGlWidget.h"
#include "MainWindow.h"

#include <QtGui/QMouseEvent>
#include <QDebug>
#include <Matrix.h>
#include <Vector.h>

MyGlWidget::MyGlWidget(QWidget *parent)
    : QGLWidget(parent), selectMode(0)
{
    setMouseTracking(true);

    timer.setInterval(200);
    timer.setSingleShot(true);
    connect(&timer, SIGNAL(timeout()), this, SLOT(mouseSingleClickEvent()));
    doubleClicked = 0;

    numPoint = 0;
    maxNumPoint = 16;
    pointIDs = new int[maxNumPoint];
    pointVertices = new GLfloat[maxNumPoint*3];
    pointColors = new GLfloat[maxNumPoint*3];

    numLine = 0;
    maxNumLine = 16;
    lineIDs = new int[maxNumLine];
    lineVertices = new GLfloat[maxNumLine*2*3];
    lineColors = new GLfloat[maxNumLine*2*3];

}

MyGlWidget::~MyGlWidget()
{

}

void MyGlWidget::setModel(MainWindow *theM)
{
    theModel = theM;
}


void
MyGlWidget::drawBuffers(){
    // draw the shapes
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);


    glLineWidth(2.0);
    glPointSize(10);

    if (numPoint > 0) {
        glColorPointer(3, GL_FLOAT, 0, pointColors);
        glVertexPointer(3, GL_FLOAT, 0, pointVertices);
        glDrawArrays(GL_POINTS, 0, numPoint);
    }

    if (numLine > 0) {
        glColorPointer(3, GL_FLOAT, 0, lineColors);
        glVertexPointer(3, GL_FLOAT, 0, lineVertices);
        glDrawArrays(GL_LINES, 0, 2*numLine);
    }
}

void
MyGlWidget::drawPoint(int tag, float x1, float y1, int numPixels, float r, float g, float b)
{
    numPoint++;
    if (numPoint > maxNumPoint) {


        GLfloat *oldPointColors = pointColors;
        GLfloat *oldPointVertices = pointVertices;
        int *oldPointIDs = pointIDs;

        int newPointSize = (maxNumPoint+32)*3;
        pointVertices = new GLfloat[newPointSize];

        pointColors = new GLfloat[newPointSize];
        pointIDs = new int[(maxNumPoint+32)];

        for (int i=0; i<maxNumPoint*3; i++) {
            pointVertices[i] = oldPointVertices[i];

            pointColors[i] = oldPointColors[i];
        }
        for (int i=0; i<maxNumPoint; i++) {
            pointIDs[i] = oldPointIDs[i];
        }

        if (oldPointVertices != 0)
            delete [] oldPointVertices;

        if (oldPointColors != 0)
            delete [] oldPointColors;
        if (oldPointIDs != 0)
            delete [] oldPointIDs;


        maxNumPoint += 32;
    }

    pointIDs[numPoint-1] = tag;

    GLfloat *locInVertices = &pointVertices[(numPoint-1)*3];
    GLfloat *locInColors = &pointColors[(numPoint-1)*3];


    // add location and value for point to pointVertices and pointColors
    locInVertices[0] = x1;
    locInVertices[1] = y1;
    locInVertices[2] = 0.;
    locInColors[0] = r;
    locInColors[1]=g;
    locInColors[2]=b;

    return;

}

void MyGlWidget::drawText(int tag, float x1, float y1, char *text, float r, float g, float b)
{
    glPushMatrix();
    glColor3f(r, g, b);
    renderText(x1, y1, 0, text);
    glPopMatrix();
}

void
MyGlWidget::drawLine(int tag, float x1, float y1, float x2, float y2, float thick, float r, float g, float b)
{
    numLine++;
    if (numLine > maxNumLine) {
        GLfloat *oldLineColors = lineColors;
        GLfloat *oldLineVertices = lineVertices;
        int *oldLineIDs = lineIDs;

        int newLineSize = (maxNumLine+32)*2*3;
        lineVertices = new GLfloat[newLineSize];
        lineColors = new GLfloat[newLineSize];
        lineIDs = new int[(maxNumLine+32)];

        for (int i=0; i<2*maxNumLine*3; i++) {
            lineVertices[i] = oldLineVertices[i];
            lineColors[i] = oldLineColors[i];
        }
        for (int i=0; i<maxNumLine; i++) {
            lineIDs[i] = oldLineIDs[i];
        }
        if (oldLineVertices != 0)
            delete [] oldLineVertices;
        if (oldLineColors != 0)
            delete [] oldLineColors;
        if (oldLineIDs != 0)
            delete [] oldLineIDs;
        maxNumLine += 32;
    }

    GLfloat *locInVertices = &lineVertices[2*(numLine-1)*3];
    GLfloat *locInColors = &lineColors[2*(numLine-1)*3];

    locInVertices[0] = x1;
    locInVertices[1] = y1;
    locInVertices[2] = 0.;
    locInVertices[3] = x2;
    locInVertices[4] = y2;
    locInVertices[5] = 0.;

    locInColors[0] = r;
    locInColors[1] = g;
    locInColors[2] = b;
    locInColors[3] = r;
    locInColors[4] = g;
    locInColors[5] = b;
}

void MyGlWidget::reset() {
    numPoint = 0;
    numLine = 0;
}

void MyGlWidget::initializeGL() {
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1, 1, 1, 0);
}

void MyGlWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//qDebug() << "HEIGHT" << buildingH;

    if (theModel != 0) {
        float heightB = theModel->getBuildingHeight();
        float maxDisp = theModel->getMaxDisp();
        float bounH = heightB/20;
        if (maxDisp == 0)
            maxDisp = 10.0;

        float bounW = maxDisp*1.1;
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
#ifdef QT_OPENGL_ES_1
        glOrthof(-maxDisp, +maxDisp, -bounH, bounH+heightB, -bounW, bounW);
#else
        glOrtho(-maxDisp, +maxDisp, -bounH, heightB+bounH, -bounW, bounW);
#endif
    } else
        glOrtho(0, 6, 0, 6, -15, 15); // set origin to bottom left corner

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void MyGlWidget::update() {

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (theModel != 0) {
        float heightB = theModel->getBuildingHeight();
        float maxDisp = theModel->getMaxDisp();
        float bounH = heightB/20;
        if (maxDisp == 0)
            maxDisp = 10.0;
        float bounW = 1.1*maxDisp;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

#ifdef QT_OPENGL_ES_1
        glOrthof(-bounW, +bounW, -bounH, bounH+heightB, -15, 15.0);
#else
        glOrtho(-bounW, +bounW, -bounH, heightB+bounH, -15, 15.0);
#endif
    } else
        glOrtho(0, 6, 0, 6, -15, 15); // set origin to bottom left corner

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    this->QGLWidget::update();
}

void MyGlWidget::paintGL() {

    glClear(GL_COLOR_BUFFER_BIT);


    if (theModel != 0)
        theModel->draw(this);

    this->drawBuffers();

}

void MyGlWidget::keyPressEvent(QKeyEvent* event) {
    switch(event->key()) {
    case Qt::Key_Escape:
        close();
        break;
    default:
        event->ignore();
        break;
    }
}

void MyGlWidget::mousePressEvent(QMouseEvent *event)
{
    // Save mouse press position
    mousePressPosition = event->localPos();

    if(event->buttons() & Qt::LeftButton) {
        clickedLeft = true;
    } else {
        clickedLeft = false;
    }
    if(doubleClicked){
        doubleClicked = 1;
    } else{
        timer.start();
    }
    QWidget::mousePressEvent(event);
}

void MyGlWidget::mouseReleaseEvent(QMouseEvent *event)
{
    mouseReleasePosition = event->localPos();

    //
    // need to determine world coords represented by mouse position in view world
    //

    GLint viewport[4]; //var to hold the viewport info
    GLdouble modelview[16]; //var to hold the modelview info
    GLdouble projection[16]; //var to hold the projection matrix info

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    for (int i=0; i<4; i++)
        viewport[i] /= devicePixelRatio();

    // get the world coordinates from the screen coordinates, in QT5.8 we have the unproject function
    // Qt5.8: QVector3D world = QVector3D(mousePressPosition.x(), mousePressPosition.y(), 0).unproject(modelview, projection, viewport, &x, &y, &z);

    Matrix viewPortMatrix(4,4);
    Matrix modelViewMatrix(modelview,4,4);
    Matrix projectionMatrix(projection, 4,4);
    Matrix A = projectionMatrix*modelViewMatrix;

    Vector in(4);
    float xPos = mousePressPosition.x();
    float yPos = viewport[3] - mousePressPosition.y() -1;

    in[0]=2*(xPos-viewport[0])/(1.0*viewport[2]) - 1.0;
    in[1]=2*(yPos-viewport[1])/(1.0*viewport[3]) - 1.0;
    in[2]=0.0;
    in[3]=1.0;
    Vector pressCrd(4);
    A.Solve(in,pressCrd);

    xPos = mouseReleasePosition.x();
    yPos = viewport[3] - mouseReleasePosition.y() -1;

    in[0]=2*(xPos-viewport[0])/(1.0*viewport[2]) - 1.0;
    in[1]=2*(yPos-viewport[1])/(1.0*viewport[3]) - 1.0;
    in[2]=0.0;
    in[3]=1.0;
    Vector releaseCrd(4);
    A.Solve(in,releaseCrd);

    // given press and release coordinated, inform ManWindow
    theModel->setSelectionBoundary(pressCrd(1),releaseCrd(1));
}

void MyGlWidget::mouseMoveEvent(QMouseEvent *event)
{
    // not used yet
}

void MyGlWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    // not used yet
    timer.stop();
    doubleClicked = 0; // this is to discard another press event coming

}

void MyGlWidget::mouseSingleClickEvent(void) {
    // not used yet
}

