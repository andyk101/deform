#include "deform_app.h"
#include "deform_dlg.h"

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <andyk/serialization/attr_xml_iarchive.hpp>
#include <andyk/serialization/attr_xml_oarchive.hpp>

#include <fstream>
#include <iomanip>
#include <locale>

// -----------------------------------------------------------------------------------------------------------
// DeformApp
// -----------------------------------------------------------------------------------------------------------
QDialog* DeformApp::findDialog(QWidget* pWidget)
{
    QDialog* pDialog = qobject_cast<QDialog*>(pWidget);
    while(pDialog == 0 && pWidget!=0)
    {
        pWidget = pWidget->parentWidget();
        pDialog = qobject_cast<QDialog*>(pWidget);
    }
    return pDialog;
}

// -----------------------------------------------------------------------------------------------------------
DeformApp::DeformApp(int& argc, char** argv)
    : QApplication(argc, argv)
{
}

DeformApp::~DeformApp()
{
}

// -----------------------------------------------------------------------------------------------------------
bool DeformApp::notify(QObject* receiver, QEvent* event)
{
    try
    {
        return QApplication::notify(receiver, event);
    }
    catch (QSharedPointer<ErrorBase> error)
    {
        QString msg = error->toString() + "\nDo you want to terminate the application?";
        QWidget* pWidget = qobject_cast<QWidget*>(receiver);
        QDialog* pDialog = findDialog(pWidget);

        int ret = QMessageBox::critical(pDialog, "Error", msg,
                                        QMessageBox::Yes|QMessageBox::Default, QMessageBox::No);
        if (ret == QMessageBox::Yes)
            exit(-1);

        return false;
    }
}

// -----------------------------------------------------------------------------------------------------------
void DeformApp::xmlParse(const XmlManager& manager)
{
    XmlParser::xmlParse(manager);
    testNames("Deform_v1.0", "Deform");
    xmlAddCollection<Material>("Material");
    xmlAddCollection<Detail>("Detail");
    xmlAddCollection<Process>("Process");
    xmlAddCollection<Criterion>("Criterion");

//    QDomElement detail = xmlChild("Detail");
//    QDomElement direction0 = xmlChild(detail, "Direction0");
//    QDomElement direction45 = xmlChild(detail, "Direction45");
//    (void)direction0;
//    (void)direction45;
}

// -----------------------------------------------------------------------------------------------------------
int DeformApp::exec()
{
    //using namespace boost::archive;
    try
    {
        Loger::instance("deform.txt");
        Loger::log("Deform was started", "", Loger::sLevel0);

//        QStringList args = arguments();
//        if (args.size() == 2 && args.at(1) == "-c")
//            validityTestDB();
//        if (args.size() == 2 && args.at(1) == "-t")
//            testTariffManager();

        // load all
//        QString source = qApp->applicationDirPath() + QDir::separator() + "deform2.xml";
//        {
//            std::ifstream ifs(source.toStdString().c_str());
//            if (!ifs.good())
//                throw ErrorBase::create(QString("Cannot open '%1' file").arg(source));
//            boost::archive::attr_xml_iarchive ia(ifs);
//            ia >> make_nvp("deform", *this);
//        }

        // config load
        XmlManager mng;
        mng.parseFile("deform.xml", this);

        // changing

        // save all
        QString target = qApp->applicationDirPath() + QDir::separator() + "deform3.xml";
        {
            std::ofstream ofs(target.toStdString().c_str());
            if (!ofs.good())
                throw ErrorBase::create(QString("Cannot open '%1' file").arg(target));
            boost::archive::attr_xml_oarchive oa(ofs);
            oa << make_nvp("deform", *this);
        }

        // деформирование
        //Process p;

        // создание gui
        DeformDlg dlg;
        dlg.resize(800, 600);

        // обработка цикла сообщений
        return dlg.exec();
    }
    catch (QSharedPointer<ErrorBase> error)
    {
        QString msg = error->toString() + "\nThe application will be terminated.";
        QMessageBox::critical(0, "Error", msg, QMessageBox::Ok);
        return -1;
    }
}

// -----------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // locale
    QLocale::setDefault(QLocale(QLocale::Russian, QLocale::RussianFederation));
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    std::locale::global(std::locale("ru_RU.UTF-8"));

#ifndef TESTING
    DeformApp app(argc, argv);
    return app.exec();
#else
    QCoreApplication app(argc, argv);
    Test_Fixture t;
    QTest::qExec(&t, app.arguments());
    return 0;
#endif
}

// -----------------------------------------------------------------------------------------------------------
