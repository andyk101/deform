#ifndef DEFORM_APP_H
#define DEFORM_APP_H

// -----------------------------------------------------------------------------------------------------------
#include <QtGui>
#include <typeinfo>
#include <andyk/core/core_fwd.h>

#include "detail.h"
//#include "deform_model.h"
//#include "deform_plot.h"
//#include "deform_table.h"

// Проводить ли тестирование
// #define TESTING

// -----------------------------------------------------------------------------------------------------------
// DeformApp
// -----------------------------------------------------------------------------------------------------------
class DeformApp;
template<class T> QMap<QString, QSharedPointer<T> >& deform_collection(DeformApp* app)
{
    throw ErrorBase::create(QString("Cannot find '%1' collection").arg(typeid(T).name()));
    static QMap<QString, QSharedPointer<T> > null;
    return null;
}

class DeformApp : public QApplication, public XmlParser {
Q_OBJECT
    // helpers
private:
    template<class T> friend QMap<QString, QSharedPointer<T> >& deform_collection(DeformApp* app);
    static QDialog* findDialog(QWidget* pWidget);

    QMap<QString, QSharedPointer<Material> >  m_materials;
    QMap<QString, QSharedPointer<Detail> >    m_details;
    QMap<QString, QSharedPointer<Process> >   m_processes;
    QMap<QString, QSharedPointer<Criterion> > m_criteria;
    QMap<QString, QSharedPointer<Curve> >     m_curves;
    QMap<QString, QSharedPointer<Plots> >     m_collsPlots;
    QMap<QString, QSharedPointer<Plot> >      m_plots;
    QMap<QString, QSharedPointer<Layouts> >   m_collsLayouts;
    QMap<QString, QSharedPointer<Layout> >    m_layouts;

//    SavingReportModel m_model;
//    SavingCostPlot m_costPlot;
//    SavingPowerPlot m_powerPlot;
//    SavingReportWeb m_reportWeb;

public:
    DeformApp(int& argc, char** argv);
    ~DeformApp();

    template<class T> QMap<QString, QSharedPointer<T> >& collection() { return deform_collection<T>(this); }
    template<class T> QSharedPointer<T> element(QString name);
    template<class T> void addCollection(QString nodeName);

    virtual void xmlParse(const XmlManager& manager);

    virtual bool notify(QObject* receiver, QEvent* event);
    virtual int exec();
};

// -----------------------------------------------------------------------------------------------------------
template<> inline QMap<QString, QSharedPointer<Material> >& deform_collection<Material>(DeformApp* app) {
    return app->m_materials;
}
template<> inline QMap<QString, QSharedPointer<Detail> >& deform_collection<Detail>(DeformApp* app) {
    return app->m_details;
}
template<> inline QMap<QString, QSharedPointer<Process> >& deform_collection<Process>(DeformApp* app) {
    return app->m_processes;
}
template<> inline QMap<QString, QSharedPointer<Criterion> >& deform_collection<Criterion>(DeformApp* app) {
    return app->m_criteria;
}
template<> inline QMap<QString, QSharedPointer<Curve> >& deform_collection<Curve>(DeformApp* app) {
    return app->m_curves;
}
template<> inline QMap<QString, QSharedPointer<Plots> >& deform_collection<Plots>(DeformApp* app) {
    return app->m_collsPlots;
}
template<> inline QMap<QString, QSharedPointer<Plot> >& deform_collection<Plot>(DeformApp* app) {
    return app->m_plots;
}
template<> inline QMap<QString, QSharedPointer<Layouts> >& deform_collection<Layouts>(DeformApp* app) {
    return app->m_collsLayouts;
}
template<> inline QMap<QString, QSharedPointer<Layout> >& deform_collection<Layout>(DeformApp* app) {
    return app->m_layouts;
}

// -----------------------------------------------------------------------------------------------------------
template<class T> inline QSharedPointer<T> DeformApp::element(QString name)
{
    QMap<QString, QSharedPointer<T> >& coll = collection<T>();
    if (!coll.contains(name))
        throw ErrorBase::create(QString("Cannot find '%1' for '%2' collection")
                                .arg(name).arg(typeid(T).name()));

    return coll[name];
}

// -----------------------------------------------------------------------------------------------------------
template<class T> void DeformApp::addCollection(QString nodeName)
{
    QMap<QString, QSharedPointer<T> >& coll = collection<T>();
    QDomElement node = xmlChildNode(xmlBody(), nodeName);
    while(!node.isNull())
    {
        QString name = xmlNodeContent(node, "Name", nodeAttr);

        bool optional = true;
        if (!coll.contains(name))
        {
            coll[name] = QSharedPointer<T>(new T());
            optional = false;
        }

        QString parent = xmlNodeContent(node, "Parent", nodeAttr, true);
        if (!parent.isEmpty())
        {
            try {
                *coll[name] = *element<T>(parent);
            } catch(QSharedPointer<ErrorBase> error) {
                throw ErrorBase::create(QString("Cannot find parent '%1' for %2 with name '%3'")
                                        .arg(parent).arg(nodeName).arg(name), error);
            }
            optional = true;
        }

        coll[name]->parse(*this, node, optional);

        node = node.nextSiblingElement(nodeName);
    }
}

// -----------------------------------------------------------------------------------------------------------
#include <QtTest>
class Test_Fixture : public QObject
{
    Q_OBJECT
public:
    Test_Fixture() {}
    virtual ~Test_Fixture() {}

private slots:
    void unit_test1()
    {
        qWarning() << "testing 1 ...";
        QVERIFY(1==1);
    }
    void unit_test2()
    {
        qWarning() << "testing 2 ...";
        QVERIFY(1==2);
    }
};

// -----------------------------------------------------------------------------------------------------------
#endif // DEFORM_APP_H
