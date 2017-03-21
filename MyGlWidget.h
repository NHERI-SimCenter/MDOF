#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QGLWidget>
#include <QVector2D>
#include <QTimer>

class MainWindow;

class MyGlWidget : public QGLWidget {

    Q_OBJECT // must include this if you use Qt signals/slots

public:
    MyGlWidget(QWidget *parent = NULL);
    ~MyGlWidget();
    void setModel(MainWindow *);
    
    void update();
    void drawLine(int tag, float x1, float y1, float x2, float y2, float thick, float r, float g, float b);
    void drawNode(int tag, float x1, float y1, int numPixels, float r, float g, float b);
    void drawText(int tag, float x1, float y1, char *text, float r, float g, float b);

  public slots:
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mouseDoubleClickEvent(QMouseEvent *);
  void mouseSingleClickEvent(void);
    
protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    /*
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    */
    void keyPressEvent(QKeyEvent *event);

    int selectMode;
    MainWindow *theModel; // pointer to class that has the draw command

    QPointF mousePressPosition;
    QPointF mouseReleasePosition;
    int doubleClicked;
    QTimer timer;
    bool clickedLeft;
};




#endif // MYGLWIDGET_H
