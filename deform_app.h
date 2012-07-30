#ifndef DEFORM_APP_H
#define DEFORM_APP_H

// -----------------------------------------------------------------------------------------------------------
#include <QtGui>
#include <typeinfo>
#include <andyk/core/core_fwd.h>
#include <andyk/serialization/serialization.h>

#include "detail.h"
//#include "deform_model.h"
//#include "deform_plot.h"
//#include "deform_table.h"

// Проводить ли тестирование
// #define TESTING

// -----------------------------------------------------------------------------------------------------------
// DeformApp
// -----------------------------------------------------------------------------------------------------------
class DeformApp : public QApplication, public XmlParser {
Q_OBJECT
    // helpers
private:
    static QDialog* findDialog(QWidget* pWidget);
    friend class boost::serialization::access;

    std::map<std::string, boost::shared_ptr<Material> > m_materials;
    std::map<std::string, boost::shared_ptr<Detail> > m_details;
    std::map<std::string, boost::shared_ptr<Process> > m_processes;
    std::map<std::string, boost::shared_ptr<Criterion> > m_criteria;

//    SavingReportModel m_model;
//    SavingCostPlot m_costPlot;
//    SavingPowerPlot m_powerPlot;
//    SavingReportWeb m_reportWeb;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
        ar & make_nvp("materials", m_materials);
        ar & make_nvp("details", m_details);
        ar & make_nvp("processes", m_processes);
        ar & make_nvp("criteria", m_criteria);
    }

public:
    DeformApp(int& argc, char** argv);
    ~DeformApp();

    std::map<std::string, boost::shared_ptr<Material> >& materials() { return m_materials; }
    std::map<std::string, boost::shared_ptr<Detail> >& details() { return m_details; }
    std::map<std::string, boost::shared_ptr<Process> >& processes() { return m_processes; }
    std::map<std::string, boost::shared_ptr<Criterion> >& criteria() { return m_criteria; }

    template<class T> std::map<std::string, boost::shared_ptr<T> >& collection();
    template<class T> boost::shared_ptr<T> element(std::string name);

    template<class T> void xmlAddCollection(QString nodeName);
    virtual void xmlParse(const XmlManager& manager);

    virtual bool notify(QObject* receiver, QEvent* event);
    virtual int exec();
};

// -----------------------------------------------------------------------------------------------------------
template<class T> inline std::map<std::string, boost::shared_ptr<T> >& deform_collection(DeformApp*) {
    throw ErrorBase::create(QString("Cannot find '%1' collection").arg(typeid(T).name()));
    static std::map<std::string, boost::shared_ptr<T> > null;
    return null;
}
template<> inline std::map<std::string, boost::shared_ptr<Material> >& deform_collection<Material>(DeformApp* app) {
    return app->materials();
}
template<> inline std::map<std::string, boost::shared_ptr<Detail> >& deform_collection<Detail>(DeformApp* app) {
    return app->details();
}
template<> inline std::map<std::string, boost::shared_ptr<Process> >& deform_collection<Process>(DeformApp* app) {
    return app->processes();
}
template<> inline std::map<std::string, boost::shared_ptr<Criterion> >& deform_collection<Criterion>(DeformApp* app) {
    return app->criteria();
}
template<class T> inline std::map<std::string, boost::shared_ptr<T> >& DeformApp::collection() {
    return deform_collection<T>(this);
}

// -----------------------------------------------------------------------------------------------------------
template<class T> inline boost::shared_ptr<T> DeformApp::element(std::string name)
{
    std::map<std::string, boost::shared_ptr<T> >& coll = collection<T>();
    if (coll.find(name) == coll.end())
        throw ErrorBase::create(QString("Cannot find '%1' for '%2' collection").arg(name.c_str()).arg(typeid(T).name()));

    return coll[name];
}

// -----------------------------------------------------------------------------------------------------------
template<class T> void DeformApp::xmlAddCollection(QString nodeName)
{
    std::map<std::string, boost::shared_ptr<T> >& coll = collection<T>();
    QDomElement node = xmlChildNode(xmlBody(), nodeName);
    while(!node.isNull())
    {
        std::string name = xmlNodeVar<std::string>(node, "Name", nodeAttr);

        bool optional = true;
        if (coll.find(name) == coll.end())
        {
            coll[name] = boost::shared_ptr<T>(new T());
            optional = false;
        }

        std::string parent = xmlNodeVar<std::string>(node, "Parent", nodeAttr, true);
        if (!parent.empty())
        {
            try {
                *coll[name] = *element<T>(parent);
            } catch(QSharedPointer<ErrorBase> error) {
                throw ErrorBase::create(QString("Cannot find parent '%1' for %2 with name '%3'")
                                        .arg(parent.c_str()).arg(nodeName).arg(name.c_str()), error);
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
