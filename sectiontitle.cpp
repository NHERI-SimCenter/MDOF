#include "sectiontitle.h"

SectionTitle::SectionTitle(QWidget *parent) : QFrame(parent)
{
    //Create Frame and Section Title
    sectionLabel = new QLabel(this);
    sectionLabel->setText(tr("Test Section"));

    // Create a section line
    line = new QFrame(); //use line instead -- needs to scale
    line->setMaximumHeight(3);
    line->setMinimumHeight(3);
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);

    //add line to Layout
    sectionLayout = new QVBoxLayout();
    sectionLayout->addWidget(sectionLabel);
    sectionLayout->addWidget(line);

    this->setLayout(sectionLayout);
}

void SectionTitle::setTitle(QString s)
{
    sectionLabel->setText(s);
}
