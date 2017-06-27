#include "SimpleSpreadsheetWidget.h"
#include <QKeyEvent>
#include <QApplication>
#include <QClipboard>
#include <QDebug>
#include <QMimeData>
#include <QMessageBox>
#include <QJsonArray>
#include <QJsonObject>

SimpleSpreadsheetWidget::SimpleSpreadsheetWidget(int colCount, int rowCount, QStringList head, QList<int>types, QWidget *parent)
    : QTableWidget(parent), numRow(rowCount), numCol(colCount), theHeadings(head), dataTypes(types)
{
    this->setRowCount(rowCount);
    this->setColumnCount(colCount);
    this->setHorizontalHeaderLabels(head);
}

SimpleSpreadsheetWidget::~SimpleSpreadsheetWidget()
{

}

void SimpleSpreadsheetWidget::clear()
{
    QTableWidget::clearContents();
}

int 
SimpleSpreadsheetWidget::getNumRows(){
  return numRow;
}
int SimpleSpreadsheetWidget::getNumColumns(){
  return numCol;
}

bool SimpleSpreadsheetWidget::getString(int row, int col, QString &res){

  QTableWidgetItem *theItem = this->item(row,col);
  if (theItem == 0)
    return false;
  res = theItem->text();
  return true;
}

bool SimpleSpreadsheetWidget::getDouble(int row, int col, double &res){

  QTableWidgetItem *theItem = this->item(row,col);
  if (theItem == 0)
    return false;
  QString textData = theItem->text();
  res = textData.toDouble();
  return true;
}

bool SimpleSpreadsheetWidget::getInt(int row, int col, int &res){
  QTableWidgetItem *theItem = this->item(row,col);
  if (theItem == 0)
    return false;
  QString textData = theItem->text();
  res = textData.toInt();
  return true;  
}

int  SimpleSpreadsheetWidget::setString(int row,int col, QString &data)
{
    QTableWidgetItem *theCell = this->item(row, col);
    if(!theCell)
    {
        theCell = new QTableWidgetItem;
        this->setItem(row, col, theCell);
    }
    theCell->setText(data);
    return 0;
}

int  SimpleSpreadsheetWidget::setDouble(int row, int col, double data)
{
    QTableWidgetItem *theCell = this->item(row, col);
    if(!theCell)
    {
        theCell = new QTableWidgetItem;
        this->setItem(row, col, theCell);
    }
    theCell->setText(QString::number(data));
    return 0;
}

int  SimpleSpreadsheetWidget::setInt(int row, int col, int data)
{
    QTableWidgetItem *theCell = this->item(row, col);
    if(!theCell)
    {
        theCell = new QTableWidgetItem;
        this->setItem(row, col, theCell);
    }
    theCell->setText(QString::number(data));
    return 0;
}

void
SimpleSpreadsheetWidget::outputToJSON(QJsonArray &rvArray)
{
   // QJsonArray rvArray;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    for (int row = 0; row < numRow; ++row) {
        QTableWidgetItem *firstItem= this->item(row, 0);
        //QModelIndex = this->it
        if (firstItem != 0) {
            QString firstItemString =firstItem->text();
            qDebug() << firstItemString;
        if (!firstItemString.isEmpty()) {
          QJsonObject obj;
          for (int column = 0; column < numCol; ++column) {
               QTableWidgetItem *theItem = this->item(row,column);
               if (theItem == 0)
                   break;
               QString textData = theItem->text();
               if (dataTypes.at(column) == SIMPLESPREADSHEET_QString)
                    obj[theHeadings.at(column)]=textData;
               else if (dataTypes.at(column) == SIMPLESPREADSHEET_QDouble) {
                   // QString textPrecision = QString("%1").arg(textData, 0, 'g', 13);
                    obj[theHeadings.at(column)]=textData.toDouble();
               } else if (dataTypes.at(column) == SIMPLESPREADSHEET_QInt)
                       obj[theHeadings.at(column)]=textData.toInt();

           }
           rvArray.append(obj);
        }

        }
    }
    QApplication::restoreOverrideCursor();
 }

void
SimpleSpreadsheetWidget::inputFromJSON(QJsonArray &rvArray){

}

void SimpleSpreadsheetWidget::keyPressEvent( QKeyEvent *event )
{

    if (event->key() == Qt::Key_C && event->modifiers() & Qt::ControlModifier) {

         QApplication::clipboard()->setText( this->currentIndex().data().toString() );
        //this->copy();
    }
    else if (event->modifiers() & Qt::ControlModifier && event->key() == Qt::Key_V) {

      model()->setData( currentIndex(), QApplication::clipboard()->text() );
        //this->paste();
    }
    else {

        QTableWidget::keyPressEvent(event);

    }
}

void SimpleSpreadsheetWidget::copy()
{
    QItemSelectionModel * selection = selectionModel();
    QModelIndexList indexes = selection->selectedIndexes();
qDebug() << "COPY";
    if(indexes.size() < 1)
        return;

    // QModelIndex::operator < sorts first by row, then by column.
    // this is what we need
//    std::sort(indexes.begin(), indexes.end());
    qSort(indexes);

    // You need a pair of indexes to find the row changes
    QModelIndex previous = indexes.first();
    indexes.removeFirst();
    QString selected_text_as_html;
    QString selected_text;
    selected_text_as_html.prepend("<html><style>br{mso-data-placement:same-cell;}</style><table><tr><td>");
    QModelIndex current;
    Q_FOREACH(current, indexes)
    {
        QVariant data = model()->data(previous);
        QString text = data.toString();
        selected_text.append(text);
        text.replace("\n","<br>");
        // At this point `text` contains the text in one cell
        selected_text_as_html.append(text);

        // If you are at the start of the row the row number of the previous index
        // isn't the same.  Text is followed by a row separator, which is a newline.
        if (current.row() != previous.row())
        {
            selected_text_as_html.append("</td></tr><tr><td>");
            selected_text.append(QLatin1Char('\n'));
        }
        // Otherwise it's the same row, so append a column separator, which is a tab.
        else
        {
            selected_text_as_html.append("</td><td>");
            selected_text.append(QLatin1Char('\t'));
        }
        previous = current;
    }

    // add last element
    selected_text_as_html.append(model()->data(current).toString());
    selected_text.append(model()->data(current).toString());
    selected_text_as_html.append("</td></tr>");
    QMimeData * md = new QMimeData;
    md->setHtml(selected_text_as_html);
//    qApp->clipboard()->setText(selected_text);
    md->setText(selected_text);
    qApp->clipboard()->setMimeData(md);
qDebug() << selected_text;
QApplication::clipboard()->setText(selected_text);
//    selected_text.append(QLatin1Char('\n'));
//    qApp->clipboard()->setText(selected_text);
}

void SimpleSpreadsheetWidget::paste()
{
    qDebug() << "PASTE";
    if(qApp->clipboard()->mimeData()->hasHtml())
    {
        // TODO, parse the html data
    }
    else
    {
        QString selected_text = qApp->clipboard()->text();
        selected_text = QApplication::clipboard()->text();
        qDebug() << selected_text;

        QStringList cells = selected_text.split(QRegExp(QLatin1String("\\n|\\t")));
        while(!cells.empty() && cells.back().size() == 0)
        {
            cells.pop_back(); // strip empty trailing tokens
        }
        int rows = selected_text.count(QLatin1Char('\n'));
        int cols = cells.size() / rows;
        if(cells.size() % rows != 0)
        {
            // error, uneven number of columns, probably bad data
            QMessageBox::critical(this, tr("Error"),
                                  tr("Invalid clipboard data, unable to perform paste operation."));
            return;
        }

        if(cols != columnCount())
        {
            // error, clipboard does not match current number of columns
            QMessageBox::critical(this, tr("Error"),
                                  tr("Invalid clipboard data, incorrect number of columns."));
            return;
        }

        // don't clear the grid, we want to keep any existing headers
        setRowCount(rows);
        // setColumnCount(cols);
        int cell = 0;
        for(int row=0; row < rows; ++row)
        {
            for(int col=0; col < cols; ++col, ++cell)
            {
                QTableWidgetItem *newItem = new QTableWidgetItem(cells[cell]);
                setItem(row, col, newItem);
            }
        }
    }
}
