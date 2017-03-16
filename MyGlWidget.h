// myglwidget.h

#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QGLWidget>
#include <QVector2D>
#include <QTimer>

class MainWindow;

class MyGlWidget : public QGLWidget
{
  Q_OBJECT
    
    public:
  
  MyGlWidget(QWidget *parent = 0);
  ~MyGlWidget();
  void setModel(MainWindow *);

  void drawLine(int tag, float x1, float y1, float x2, float y2, float thick, float r, float g, float b);
  void drawNode(int tag, float x1, float y1, int numPixels, float r, float g, float b);

 protected:
  void initializeGL();
  void paintGL();
  void resizeGL(int width, int height);
  
  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  public slots:
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *);
  void mouseSingleClickEvent(void);
  
 private:
  void draw();
  QPoint lastPos;

  GLfloat *verts;
  GLuint vertexBufferID;

  GLshort *indices;
  GLuint eleBufferID;

  MainWindow *theModel; // pointer to class that has the draw command
  int selectMode;

  QPointF mousePressPosition;
  QPointF mouseReleasePosition;
  int doubleClicked;
  QTimer timer;
  bool clickedLeft;
};

#endif // MYGLWIDGET_H
