#ifndef SIMPLESPREADSHEETWIDGET_H
#define SIMPLESPREADSHEETWIDGET_H

#include <QTableWidget>
#include <QJsonArray>
#include <QStringList>
#include <QList>
#include <QString>

#define SIMPLESPREADSHEET_QString 0
#define SIMPLESPREADSHEET_QDouble 1
#define SIMPLESPREADSHEET_QInt    2

class QKeyEvent;


class SimpleSpreadsheetWidget : public QTableWidget
{
    Q_OBJECT
public:
    explicit SimpleSpreadsheetWidget(int colCount, int rowCount, QStringList, QList<int>, QWidget *parent = 0);
    ~SimpleSpreadsheetWidget();

    void outputToJSON(QJsonArray &rvArray);
    void inputFromJSON(QJsonArray &rvArray);

    int getNumRows();
    int getNumColumns();
    bool getString(int row, int col, QString &);
    bool getDouble(int row, int col, double &);
    bool getInt(int row, int col, int &);
    int  setString(int row,int col, QString &);
    int  setDouble(int row, int col, double);
    int  setInt(int row, int col, int);
    
protected:
    void keyPressEvent( QKeyEvent *event );

signals:

public slots:
    void clear(void);

private:
    void copy();
    void paste();
    int numRow;
    int numCol;

    QStringList theHeadings;
    QList<int>  dataTypes;
};

#endif // SIMPLESPREADSHEETWIDGET_H
