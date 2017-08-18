#ifndef SECTIONTITLE_H
#define SECTIONTITLE_H

#include <QWidget>
#include <QString>
#include <QFrame>
#include <QLabel>
#include <QVBoxLayout>

class SectionTitle : public QFrame
{
    Q_OBJECT
public:
    explicit SectionTitle(QWidget *parent = 0);
    void setTitle(QString);

signals:

public slots:

private:
    QVBoxLayout *sectionLayout;
    QLabel      *sectionLabel;
    QFrame      *line;
};

#endif // SECTIONTITLE_H
