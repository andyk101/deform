#include "deform_app.h"
#include "deform_dlg.h"

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
    addCollection<Material>("Material");
    addCollection<Detail>("Detail");
    addCollection<Process>("Process");
    addCollection<Criterion>("Criterion");
    addCollection<Curve>("Curve");
    addCollection<Plots>("Plots");
    addCollection<Layouts>("Layouts");

    foreach(QSharedPointer<Plots> plots, m_collsPlots.values())
    {
        plots->create();
    }
    foreach(QSharedPointer<Layouts> layouts, m_collsLayouts.values())
    {
        layouts->create();
    }
}

// -----------------------------------------------------------------------------------------------------------
int DeformApp::exec()
{
    //using namespace boost::archive;
    //try
    {
        Loger::instance("deform.txt");
        Loger::log("Deform was started", "", Loger::sLevel0);

//        QStringList args = arguments();
//        if (args.size() == 2 && args.at(1) == "-c")
//            validityTestDB();
//        if (args.size() == 2 && args.at(1) == "-t")
//            testTariffManager();

        // config load
        XmlManager mng;
        mng.parseFile("deform.xml", this);

        // деформирование
        //Plot::find("m_d(r_km)-08кп-Local")->calculate();
        QSharedPointer<Process> process = Plot::find("m_d(r_km)-08кп-Local")->criterion()->process();
        process->exec();

        // создание gui
        DeformDlg dlg;
        dlg.resize(800, 600);

        // обработка цикла сообщений
        return dlg.exec();
    }
//    catch (QSharedPointer<ErrorBase> error)
//    {
//        QString msg = error->toString() + "\nThe application will be terminated.";
//        QMessageBox::critical(0, "Error", msg, QMessageBox::Ok);
//        return -1;
//    }
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
