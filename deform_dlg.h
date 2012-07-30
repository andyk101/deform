#ifndef DEFORM_DLG_H
#define DEFORM_DLG_H

// -----------------------------------------------------------------------------------------------------------
#include <QtGui>

// -----------------------------------------------------------------------------------------------------------
// DeformDlg
// -----------------------------------------------------------------------------------------------------------
class DeformDlg : public QDialog {
    Q_OBJECT
public:
    // static

private:
    QGridLayout* m_pLayout;
    QPushButton* m_pbtnOk;

public:
    DeformDlg(QWidget* parent = 0);

protected:
    virtual void resizeEvent ( QResizeEvent * );

public slots:
    void ok();
};

// -----------------------------------------------------------------------------------------------------------
#endif // DEFORM_DLG_H
