#include "sectiontitle.h"

SectionTitle::SectionTitle(QWidget *parent) : QFrame(parent)
{
    //
    //
    //Create Frame and Section Title


    sectionLayout = new QVBoxLayout();
    this->setLayout(sectionLayout);
    sectionLabel = new QLabel(this);
    titleText->setText(tr("Test Section"));
    //add title to Layout


    //
    // Create a section line

    line = new QFrame(); //use line instead -- needs to scale
    line->setGeometry(QRect(320, 150, 118, 3)); //set min and max height
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    //add line to Layout


}

//sectionTitle::setText
//where you can change the text
//void function
