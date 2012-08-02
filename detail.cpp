#include "detail.h"
#include "deform_app.h"

// -----------------------------------------------------------------------------------------------------------
// Material - материал
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Material> Material::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Material>(name);
}

void Material::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name,element, "Name",XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_R0,   element, "R0",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_R45,  element, "R45", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_B,    element, "B",   XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_m,    element, "m",   XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_Omega,element, "Omega", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_U,    element, "U",   XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_aa0,  element, "aa0", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_aa1,  element, "aa1", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_aa2,  element, "aa2", XmlParser::nodeAttr, optional);
}

// -----------------------------------------------------------------------------------------------------------
// Detail
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Detail> Detail::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Detail>(name);
}

void Detail::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    QString materialName = parser.xmlNodeVar<QString>(element, "Material",  XmlParser::nodeAttr, optional);
    if (!materialName.isEmpty())
        m_material = Material::find(materialName);

    parser.xmlNodeVarChange<double>(m_r_0,  element, "r_0",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_m_d,  element, "m_d",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_r_kp, element, "r_kp", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_r_km, element, "r_km", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_s_0,  element, "s_0",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_z,    element, "z",    XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_mu,   element, "mu",   XmlParser::nodeAttr, optional);
}

QSharedPointer<Detail> Detail::clone() const
{
    QSharedPointer<Detail> det(new Detail(*this));
    det->m_name = "";
    return det;
}

// -----------------------------------------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Process> Process::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Process>(name);
}

Process& Process::operator= (const Process& proc)
{
    m_name = proc.m_name;
    m_detail = proc.m_detail->clone();
    m_v_parts = proc.m_v_parts;
    m_calc_s = proc.m_calc_s;
    return *this;
}

void Process::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    QString detailName = parser.xmlNodeVar<QString>(element, "Detail",  XmlParser::nodeAttr, optional);
    if (!detailName.isEmpty())
        m_detail = Detail::find(detailName)->clone();

    parser.xmlNodeVarChange<int>(m_v_parts, element, "v_parts", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<bool>(m_calc_s, element, "calc_s",  XmlParser::nodeAttr, optional);
}

QSharedPointer<Process> Process::clone() const
{
    QSharedPointer<Process> proc(new Process(*this));
    proc->m_name = "";
    return proc;
}

// -----------------------------------------------------------------------------------------------------------
// Critarion
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Criterion> Criterion::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Criterion>(name);
}

Criterion& Criterion::operator= (const Criterion& crit)
{
    m_name = crit.m_name;
    m_process = crit.m_process->clone();
    m_type = crit.m_type;
    m_method = crit.m_method;
    return *this;
}

void Criterion::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    QString processName = parser.xmlNodeVar<QString>(element, "Process",  XmlParser::nodeAttr, optional);
    if (!processName.isEmpty())
        m_process = Process::find(processName)->clone();

    parser.xmlNodeVarChange<CriterionType>(m_type, element, "Type", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_method, element, "Method", XmlParser::nodeAttr, optional);
}

QSharedPointer<Criterion> Criterion::clone() const
{
    QSharedPointer<Criterion> crit(new Criterion(*this));
    crit->m_name = "";
    return crit;
}

// -----------------------------------------------------------------------------------------------------------
// Curve
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Curve> Curve::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Curve>(name);
}

void Curve::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name, element, "Name", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<CurveArgType>(m_argType, element, "ArgType", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_argStart, element, "ArgStart", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_argEnd, element, "ArgEnd", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_argStep, element, "ArgStep", XmlParser::nodeAttr, optional);
}

// -----------------------------------------------------------------------------------------------------------
// Plots
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Plots> Plots::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Plots>(name);
}

QString Plots::cartesianName(QString strTemplate, QStringList strCartesian)
{
    QString name = strTemplate;
    for (PlotsCollection coll(0); coll.isValid(); ++coll)
    {
        QString pattern = QString("{%1}").arg((QString)coll);
        name.replace(pattern, strCartesian[(int)coll]);
    }
    return name;
}

void Plots::addString(QDomElement& elem, QStringList* lst)
{
    QString item = elem.text();

    if (!lst->contains(item))
        *lst << item;
}

void Plots::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name,     element, "Name",     XmlParser::nodeAttr,  optional);
    parser.xmlNodeVarChange<QString>(m_plotName, element, "PlotName", XmlParser::nodeChild, optional);

    for (PlotsCollection coll(0); coll.isValid(); ++coll)
    {
        bool currOptional = ((int)coll == PlotsCollections::eMaterial ||
                             (int)coll == PlotsCollections::eDetail ||
                             (int)coll == PlotsCollections::eProcess) ? true : optional;
        QStringList& lst = m_collections[(int)coll];
        parser.xmlParseCollection(this, &Plots::addString, &lst, element, coll, currOptional);
    }
}

void Plots::create()
{
    QMap<QString, QSharedPointer<Plot> >& plots = qobject_cast<DeformApp*>(qApp)->collection<Plot>();

    Cartesian cartesian(PlotsCollections::eCount);
    for (PlotsCollection coll(0); coll.isValid(); ++coll)
    {
        int count = m_collections[(int)coll].count();
        cartesian.setLimit((int)coll, (count > 1) ? count : 1);
    }

    for(cartesian.first(); cartesian.isValid(); cartesian.next())
    {
        QStringList strCartesian;
        for (PlotsCollection coll(0); coll.isValid(); ++coll)
        {
            QStringList& lst = m_collections[(int)coll];
            strCartesian << (lst.isEmpty() ? QString() : lst[ cartesian[(int)coll] ]);
        }

        QString name = cartesianName(m_plotName, strCartesian);
        if (!plots.contains(name))
        {
            plots[name] = QSharedPointer<Plot>(new Plot());
        }

        plots[name]->reset(name, strCartesian);
    }
}

// -----------------------------------------------------------------------------------------------------------
// Plot
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Plot> Plot::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Plot>(name);
}

void Plot::reset(QString name, QStringList strCartesian)
{
    m_name = name;
    m_criterion = Criterion::find(strCartesian[PlotsCollections::eCriterion])->clone();
    m_curve = Curve::find(strCartesian[PlotsCollections::eCurve]);
    m_empty = true;

    QString strProc = strCartesian[PlotsCollections::eProcess];
    if (!strProc.isEmpty())
        m_criterion->setProcess(Process::find(strProc)->clone());

    QString strDetail = strCartesian[PlotsCollections::eDetail];
    if (!strDetail.isEmpty())
        m_criterion->process()->setDetail(Detail::find(strDetail)->clone());

    QString strMaterial = strCartesian[PlotsCollections::eMaterial];
    if (!strMaterial.isEmpty())
        m_criterion->process()->detail()->setMaterial(Material::find(strMaterial));
}

void Plot::calculate()
{
    // ...
    m_empty = false;
}

// -----------------------------------------------------------------------------------------------------------
// Layouts
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Layouts> Layouts::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Layouts>(name);
}

void Layouts::addString(QDomElement& elem, QStringList* lst)
{
    QString item = elem.text();

    if (!lst->contains(item))
        *lst << item;
}

void Layouts::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<QString>(m_name,     element, "Name",     XmlParser::nodeAttr,  optional);
    parser.xmlNodeVarChange<QString>(m_layoutName, element, "LayoutName", XmlParser::nodeChild, optional);

    parser.xmlParseCollection(this, &Layouts::addString, &m_plots, element, "Plot", optional);
    for (PlotsCollection coll(0); coll.isValid(); ++coll)
    {
        QStringList* lst = &m_collections[(int)coll];
        parser.xmlParseCollection(this, &Layouts::addString, lst, element, coll, true);
    }
}

void Layouts::create()
{
    QMap<QString, QSharedPointer<Layout> >& layouts = qobject_cast<DeformApp*>(qApp)->collection<Layout>();

    Cartesian cartesian(PlotsCollections::eCount);
    for (PlotsCollection coll(0); coll.isValid(); ++coll)
    {
        int count = m_collections[(int)coll].count();
        cartesian.setLimit((int)coll, (count > 1) ? count : 1);
    }

    for(cartesian.first(); cartesian.isValid(); cartesian.next())
    {
        QStringList strCartesian;
        for (PlotsCollection coll(0); coll.isValid(); ++coll)
        {
            QStringList& lst = m_collections[(int)coll];
            strCartesian << (lst.isEmpty() ? QString() : lst[ cartesian[(int)coll] ]);
        }

        QString name = Plots::cartesianName(m_layoutName, strCartesian);
        if (!layouts.contains(name))
        {
            layouts[name] = QSharedPointer<Layout>(new Layout());
        }

        QStringList plotNames;
        foreach(QString plot, m_plots)
        {
            plotNames << Plots::cartesianName(plot, strCartesian);
        }

        layouts[name]->reset(name, plotNames);
    }
}

// -----------------------------------------------------------------------------------------------------------
// Layout
// -----------------------------------------------------------------------------------------------------------
QSharedPointer<Layout> Layout::find(QString name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Layout>(name);
}

void Layout::reset(QString name, QStringList plotNames)
{
    m_name = name;
    m_plots.clear();
    foreach(QString plotName, plotNames)
    {
        m_plots << Plot::find(plotName);
    }
    m_empty = true;
}

void Layout::create()
{
    // ...
    m_empty = false;
}


// -----------------------------------------------------------------------------------------------------------
