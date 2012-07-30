#include "detail.h"
#include "deform_app.h"

// -----------------------------------------------------------------------------------------------------------
// Material - материал
// -----------------------------------------------------------------------------------------------------------
boost::shared_ptr<Material> Material::find(std::string name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Material>(name);
}

void Material::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<std::string>(m_name, element, "Name", XmlParser::nodeAttr, optional);
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
boost::shared_ptr<Detail> Detail::find(std::string name)
{
    return qobject_cast<DeformApp*>(qApp)->element<Detail>(name);
}

void Detail::parse(XmlParser& parser, QDomElement& element, bool optional)
{
    parser.xmlNodeVarChange<std::string>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    std::string materialName = parser.xmlNodeVar<std::string>(element, "Material",  XmlParser::nodeAttr, optional);
    if (!materialName.empty())
        m_material = Material::find(materialName);

    parser.xmlNodeVarChange<double>(m_r_0,  element, "r_0",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_m_d,  element, "m_d",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_r_kp, element, "r_kp", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_r_km, element, "r_km", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_s_0,  element, "s_0",  XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_z,    element, "z",    XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_mu,   element, "mu",   XmlParser::nodeAttr, optional);
}

boost::shared_ptr<Detail> Detail::clone() const
{
    boost::shared_ptr<Detail> det(new Detail(*this));
    det->m_name = "";
    return det;
}

// -----------------------------------------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------------------------------------
boost::shared_ptr<Process> Process::find(std::string name)
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
    parser.xmlNodeVarChange<std::string>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    std::string detailName = parser.xmlNodeVar<std::string>(element, "Detail",  XmlParser::nodeAttr, optional);
    if (!detailName.empty())
        m_detail = Detail::find(detailName)->clone();

    parser.xmlNodeVarChange<int>(m_v_parts, element, "v_parts", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<bool>(m_calc_s, element, "calc_s",  XmlParser::nodeAttr, optional);
}

boost::shared_ptr<Process> Process::clone() const
{
    boost::shared_ptr<Process> proc(new Process(*this));
    proc->m_name = "";
    return proc;
}

// -----------------------------------------------------------------------------------------------------------
// Critarion
// -----------------------------------------------------------------------------------------------------------
boost::shared_ptr<Criterion> Criterion::find(std::string name)
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
    parser.xmlNodeVarChange<std::string>(m_name, element, "Name", XmlParser::nodeAttr, optional);

    std::string processName = parser.xmlNodeVar<std::string>(element, "Process",  XmlParser::nodeAttr, optional);
    if (!processName.empty())
        m_process = Process::find(processName)->clone();

    parser.xmlNodeVarChange<CriterionType>(m_type, element, "Type", XmlParser::nodeAttr, optional);
    parser.xmlNodeVarChange<double>(m_method, element, "Method", XmlParser::nodeAttr, optional);

    m_type = CriterionTypes::eFenom;
    if (m_type == CriterionTypes::eFenom)
    {
        qWarning() << (QString)m_type;
    }
    m_type = "Fenom";
    if ((QString)m_type == "Fenom")
    {
        qWarning() << (int)m_type;
    }

    int i = sizeof(m_type);
}

boost::shared_ptr<Criterion> Criterion::clone() const
{
    boost::shared_ptr<Criterion> crit(new Criterion(*this));
    crit->m_name = "";
    return crit;
}

// -----------------------------------------------------------------------------------------------------------
// Plot - график
// -----------------------------------------------------------------------------------------------------------
Plot::Plot(QWeakPointer<Detail> detal, QWeakPointer<Detail> curve, QWeakPointer<Detail> material,
     QWeakPointer<Criterion> critarion)
{
    (void)detal;
    (void)curve;
    (void)material;
    (void)critarion;
}

// -----------------------------------------------------------------------------------------------------------
