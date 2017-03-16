#include <QtWidgets>
#include <QtOpenGL>

#include "myglwidget.h"
#include "MainWindow.h"
#include <QVector3D>
#include <Matrix.h>
#include <Vector.h>

MyGlWidget::MyGlWidget(QWidget *parent)
  : QGLWidget(QGLFormat(QGL::SampleBuffers), parent), selectMode(0)
{
  timer.setInterval(200);
  timer.setSingleShot(true);
  connect(&timer, SIGNAL(timeout()), this, SLOT(mouseSingleClickEvent()));
  doubleClicked = 0;
}

MyGlWidget::~MyGlWidget()
{
}

void MyGlWidget::setModel(MainWindow *theM)
{
    theModel = theM;
}

QSize MyGlWidget::minimumSizeHint() const
{
  return QSize(50, 50);
}

QSize MyGlWidget::sizeHint() const
{
  return QSize(400, 400);
}


void MyGlWidget::initializeGL()
{
  qglClearColor(Qt::white);

  glEnable(GL_MULTISAMPLE);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glShadeModel(GL_SMOOTH);
  //glEnable(GL_LIGHTING);
  //glEnable(GL_LIGHT0);

  //static GLfloat lightPosition[4] = { 0, 0, 10, 1.0 };
  //glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

void MyGlWidget::paintGL()
{
  double xRot = 0;
  double yRot = 0;
  double zRot = 0;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glLoadIdentity();
 // glTranslatef(0.0, 0.0, -10.0);
 // glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
 // glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
 // glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);

  if (theModel != 0) {
      float heightB = theModel->getHeight();
      //heightB = 10;
      float maxDisp = theModel->getMaxDisp();
      float bounH = heightB/20;
      if (maxDisp == 0)
              maxDisp = 10.0;
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
    #ifdef QT_OPENGL_ES_1
      glOrthof(-maxDisp, +maxDisp, -bounH, bounH+heightB, -15, 15.0);
    #else
       glOrtho(-maxDisp, +maxDisp, -bounH, heightB+bounH, -15, 15.0);
    #endif
//qDebug() << maxDisp << " " << bounH << " " << heightB;
  }

  glMatrixMode(GL_MODELVIEW);

  draw();
}

void MyGlWidget::resizeGL(int width, int height)
{
   height = this->height();
   width = this->width();
  int side = qMin(width, height);
 // glViewport((width - side) / 2, (height - side) / 2, side, side);

  glViewport(0, 0, side, side);
  qDebug() << "side: " << side << " " << width << " " << height;

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

    float xPos = mouseReleasePosition.x();
    float yPos = viewport[3] - mouseReleasePosition.y() -1;

    in[0]=2*(xPos-viewport[0])/(1.0*viewport[2]) - 1.0;
    in[1]=2*(yPos-viewport[1])/(1.0*viewport[3]) - 1.0;
    in[2]=0.0;
    in[3]=1.0;
    Vector releaseCrd(4);
    A.Solve(in,releaseCrd);

    // given press and release coordinated, set vectors containing floors and stories that can be edited
    theModel->setSelectionBoundary()

}

void MyGlWidget::mouseMoveEvent(QMouseEvent *event)
{
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();

  if (event->buttons() & Qt::LeftButton) {
    //    setXRotation(xRot + 8 * dy);
    //    setYRotation(yRot + 8 * dx);
  } else if (event->buttons() & Qt::RightButton) {
    //    setXRotation(xRot + 8 * dy);
    //    setZRotation(zRot + 8 * dx);
  }

  lastPos = event->pos();
}

void MyGlWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
  qDebug() << clickedLeft << " DOUBE CLICK" << doubleClicked;

  timer.stop();
  doubleClicked = 0; // this is to discard another press event coming
  
  lastPos = event->pos();


}

void MyGlWidget::mouseSingleClickEvent(void) {
 qDebug() << clickedLeft << " SINGLE CLICK ";
}

void
MyGlWidget::drawNode(int tag, float x1, float y1, int numPixels, float r, float g, float b)
{
   // qDebug() << x1 << " " << y1 << " " << numPixels << r;
    if (selectMode == 0) {
        glPointSize(numPixels);
        glColor3f(r, g, b);
        glBegin(GL_POINTS);
        glVertex3f(x1, y1, 0.0);
        glEnd();
    } else {
        if (tag != 0) {
            glPointSize(numPixels);
            int r1 = (tag & 0x000000FF) >>  0;
            int g1 = (tag & 0x0000FF00) >>  8;
            int b1 = (tag & 0x00FF0000) >> 16;
                   glColor3f(r1/255.0,g1/255.0,b1/255.0);
    glPointSize(numPixels);
        glBegin(GL_POINTS);
        glVertex3f(x1, y1, 0.0);
        glEnd();
        }
    }
}

void
MyGlWidget::drawLine(int tag, float x1, float y1, float x2, float y2, float thick, float r, float g, float b)
{
    if (selectMode == 0) {
        glLineWidth(thick);
        glColor3f(r, g, b);
        glBegin(GL_LINES);
        glVertex3f(x1, y1, 0.0);
        glVertex3f(x2, y2, 0.0);
        glEnd();
    } else {
        if (tag != 0) {
            int r1 = (tag & 0x000000FF) >>  0;
            int g1 = (tag & 0x0000FF00) >>  8;
            int b1 = (tag & 0x00FF0000) >> 16;
                   glColor3f(r1/255.0,g1/255.0,b1/255.0);
        glLineWidth(thick);
        glBegin(GL_LINES);
        glVertex3f(x1, y1, 0.0);
        glVertex3f(x2, y2, 0.0);
        glEnd();
        }
    }
}

void MyGlWidget::draw()
{
    if (theModel != 0)
        return theModel->draw(this);
}
