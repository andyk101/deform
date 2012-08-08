#include "detail.h"
#include "deform_app.h"

// ----------------------------------------------------------------------------------------------------------
// Решение квадратного уравнения a*x^2+b*x+c=0, если a,b,c - числа с плав. точкой.
// ----------------------------------------------------------------------------------------------------------
bool SquareSolve(double& x1, double& x2, double a, double b, double c)
{
    if (a == 0)
    {
        // c=0
        if (b == 0)
            return false;

        // b*x+c=0
        x1 = x2 = -c/b;
        return true;
    }

    // x^2+b*x+c=0
    b /= a;
    c /= a;

    double d = b*b-4*c;
    if (d < 0)
        return false;

    d = sqrt(d)/2;
    b = -b/2;
    x1 = b-d;
    x2 = b+d;
    return true;
}

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

void Detail::first_h(int v_parts)
{
    m_r_2 = m_m_d*m_r_0 - m_z/2;
    m_r_c = m_r_km + m_z + m_r_2;
    m_r_max = m_m_d*m_r_0;
    m_V0 = M_PI*m_s_0*pow(m_r_0, 2);
    m_V1 = M_PI*m_s_0*pow(m_r_2-m_r_kp, 2);
    m_v_parts = v_parts;
//    double V7 = M_PI*m_s_0*( pow(m_r_0,2) - pow(m_r_c,2) );
//    m_v_parts = (int)round(V7 / m_V0 * v_parts);

    m_geom[eDeg0] = QSharedPointer<Geom>(new Geom(this, eDeg0));
    m_geom[eDeg45] = QSharedPointer<Geom>(new Geom(this, eDeg45));
}

bool Detail::isValid() const
{
    return m_geom[eDeg0]->isValid() && m_geom[eDeg45]->isValid();
}

void Detail::next_h(double dh)
{
    // Возможно следует поднять m_h выше???
    m_geom[eDeg0]->next_h(dh);
    m_geom[eDeg45]->next_h(dh);
}

// -----------------------------------------------------------------------------------------------------------
// Point
// -----------------------------------------------------------------------------------------------------------
GeomsPoint::GeomsPoint()
{
    throw ErrorBase::create(QString("GeomsPoint() is invalid constructor"));
}

GeomsPoint::GeomsPoint(Geom* geom, double s)
    : m_geom(geom), m_s(s)
{
    m_v = m_r = m_h = m_epsilon_phi = m_epsilon_i =
          m_sigma_s = m_sigma_r = m_sigma_phi = m_s_expr = m_omega_e = 0;
}

// -----------------------------------------------------------------------------------------------------------
// Geom
// -----------------------------------------------------------------------------------------------------------
Geom::Geom(Detail* detail, int direction)
    : m_material(&*detail->material()), m_detail(detail), m_direction(direction),
      m_points(v_parts()+1, GeomsPoint(this, s_0())), m_bounds(eBoundCount, GeomsPoint(this, s_0())), m_valid(false)
{
    m_h = 0;
    m_s_1 = s_0();
    m_V2 = m_V6 = 0;
    m_V5 = M_PI*s_0()*( pow(r_c(),2) - pow(r_2()-r_kp(),2) );
    m_V7 = M_PI*s_0()*( pow(r_0(),2) - pow(r_c(),2) );
    m_alpha = 0;
    m_AB = r_c() - (r_2() - r_kp());
    m_r_1 = r_c();
    m_r_k = r_0();

    m_V7_i_max = m_V6_i_max = m_V5_i_max = -1;
    for (int i = 0; i < m_points.count(); i++)
    {
        GeomsPoint& pt = m_points[i];
        pt.m_v = V0() - i*dv();

        pt.m_r = (i == 0) ? r_0() :
          (i < v_parts()) ? sqrt(pow(m_points[i-1].m_r,2) - dv()/(M_PI*s_0())) :
                            0;

        if (pt.m_r >= r_c())
        {
            m_V7_i_max = m_V6_i_max = m_V5_i_max = i;
        }
        else if (pt.m_r >= r_2() - r_kp())
        {
            m_V5_i_max = i;
        }
        else
        {
            break;
        }
    }

    m_bounds[eBoundV6].m_v = m_bounds[eBoundV7].m_v = V0() - m_V7;
    m_bounds[eBoundV6].m_r = m_bounds[eBoundV7].m_r = r_c();
    m_bounds[eBoundRmax].m_v = M_PI*s_0()*pow(r_max(),2);
    m_bounds[eBoundRmax].m_r = r_max();
    m_bounds[eBoundV5].m_v = V1();
    m_bounds[eBoundV5].m_r = r_2() - r_kp();
    m_valid = true;
}

// -----------------------------------------------------------------------------------------------------------
double Geom::s_1_x(double AB_x) const
{
    if (AB_x < 0 || AB_x > m_AB )
        throw ErrorBase::create(QString("AB_x must be in range [0, %1], but actual is %2").arg(m_AB).arg(AB_x));
    return s_0()+(m_s_1-s_0())*AB_x/m_AB;
}

double Geom::V2_x(double alpha_x) const
{
    if (alpha_x < 0 || alpha_x > m_alpha )
        throw ErrorBase::create(QString("alpha_x must be in range [0, %1], but actual is %2").arg(m_alpha).arg(alpha_x));
    return 4*M_PI/3*s_0()*pow(sin(alpha_x/2),2)*(3*pow(r_kp(),2) + 3*r_kp()*s_0() + pow(s_0(),2)) +
        M_PI*alpha_x*s_0()*(2*r_kp()+s_0())*(r_2()-r_kp());
}

double Geom::V5_x(double AB_x) const
{
    if (AB_x < 0 || AB_x > m_AB )
        throw ErrorBase::create(QString("AB_x must be in range [0, %1], but actual is %2").arg(m_AB).arg(AB_x));
    double s_1_xx = s_1_x(AB_x);
    return M_PI*AB_x/(3*cos(m_alpha)) * ( AB_x*(s_0()+2*s_1_xx) + 3*(s_0()+s_1_xx)*(r_2()-r_kp()) +
        sin(m_alpha)*((s_0()+s_1_xx)*(3*r_kp()+2*s_0())-pow(s_1_xx,2)) );
}

double Geom::V6_x(double alpha_x) const
{
    if (alpha_x < 0 || alpha_x > m_alpha )
        throw ErrorBase::create(QString("alpha_x must be in range [0, %1], but actual is %2").arg(m_alpha).arg(alpha_x));
    return M_PI*alpha_x*m_s_1*(2*r_km()+m_s_1)*r_c() - 4.0/3*M_PI*m_s_1*pow(sin(alpha_x/2),2)*
        (3*pow(r_km(),2)+3*r_km()*m_s_1+pow(m_s_1,2));
}

double Geom::V7_x(double r_k_x) const
{
    if (r_k_x < r_c() || r_k_x > m_r_k )
        throw ErrorBase::create(QString("r_k_x must be in range [%1, %2], but actual is %3").arg(r_c()).arg(m_r_k).arg(r_k_x));
    return M_PI*m_s_1*(pow(r_k_x,2) - pow(r_c(),2));
}

// -----------------------------------------------------------------------------------------------------------
void Geom::calcPoint(double& r_x, double& h_x, double v_x) const
{
    if (v_x < V1()+m_V2 || v_x > V0() )
        throw ErrorBase::create(QString("v_x must be in range [%1, %2], but actual is %3").arg(V1()+m_V2).arg(V0()).arg(v_x));

    // V7
    if (m_V7 > 0 && v_x >= V0()-m_V7)
    {
        r_x = sqrt((v_x-V0()+m_V7)/(M_PI*m_s_1+pow(r_c(),2)));
        h_x = 0;
    }
    // V6
    else if (v_x >= V1()+m_V2+m_V5)
    {
        //alpha_x:=fsolve(fV6(alpha_xx)=V1()+m_V2+m_V5+V6_x(m_alpha)-v_x, alpha_xx=0..alpha);
        //r_x = r_c()-(r_km()+s_0()/2)*sin(alpha_x);
        //h_z = (r_km()+s_0()/2)*(1-cos(alpha_x));
    }
    // V5
    else if (v_x >= V1()+m_V2)
    {
    }
}

void Geom::next_h(double dh)
{
    if (!isValid())
        return;

    m_h = m_h + dh;

    // -------------------------------------------------------------------------------------------------------
    // (m_h) -> (m_alpha, m_AB, m_r_1)

    // fh(alpha_xx)=h -> sin(alpha_x)*a1+cos(alpha_x)*b1-c1=0
    double a1 = r_c()-r_2()+r_kp();
    double b1 = r_km()+r_kp()+s_0()-m_h;
    double c1 = r_km()+r_kp()+s_0();

    if (b1 == 0.0)
    {
        m_valid = a1 != 0;
        if (!m_valid)
            return;
        m_alpha = asin(c1/a1);
    }
    else
    {
        // { sin(alpha_x)*a1+cos(alpha_x)*b1-c1=0, sin(alpha_x)^2+cos(alpha_x)^2=1 }
        // x=sin(alpha_x) -> y=cos(alpha_x)=(c1-x*a1)/b1
        // x^2+((c1-x*a1)/b1)^2=1 -> { a2*x^2+b2*x+c2=0, y=(c1-x*a1)/b, x>=0, y>=0, alpha_x=arcsin(x) }
        double a2 = pow(a1,2)+pow(b1,2);
        double b2 = -2*a1*c1;
        double c2 = pow(c1,2)-pow(b1,2);

        double x1, x2;
        m_valid = SquareSolve(x1,x2, a2,b2,c2);
        if (!m_valid)
            return;

        double x;
        if (x1>=0 && (c1-x1*a1)/b1>=0)
            x = x1;
        else if (x2>=0 && (c1-x2*a1)/b1>=0)
            x = x2;
        else
        {
            m_valid = false;
            return;
        }
        m_alpha = asin(x);
    }
    m_AB = r_c() - r_2() + r_kp() - sin(m_alpha)*(r_km()+r_kp()+s_0());
    m_r_1 = r_c() - (r_km()+s_0()/2)*sin(m_alpha);

    // -------------------------------------------------------------------------------------------------------
    // (s_0, s_1, alpha, AB, r_1) -> (V2, V5, V6, V7)
    m_V2 = V2_x(m_alpha);
    m_V5 = V5_x(m_AB);
    m_V6 = V6_x(m_alpha);
    if (V0() < V1()+m_V2+m_V5)
    {
        m_valid = false;
        return;
    }
    else if (V0() < V1()+m_V2+m_V5+m_V6)
    {
        m_V6 = V0() - (V1()+m_V2+m_V5);
        m_V7 = 0;
    }
    else
    {
        m_V7 = V0() - (V1()+m_V2+m_V5+m_V6);
    }

    // -------------------------------------------------------------------------------------------------------
    // (points_prev) -> (r_k, points, V7_i_max, V6_i_max, V5_i_max)
    int V5_i_max_prev = m_V5_i_max;
    m_V7_i_max = m_V6_i_max = m_V5_i_max = -1;
    for (int i = 0; i < V5_i_max_prev; i++)
    {
        GeomsPoint& pt = m_points[i];
        calcPoint(pt.m_r, pt.m_h, pt.m_v);

        if (i == 0)
            m_r_k = pt.m_r;

        if (pt.m_v >= V0() - m_V7)
        {
            m_V7_i_max = m_V6_i_max = m_V5_i_max = i;
        }
        else if (pt.m_v >= V0() - (m_V7+m_V6))
        {
            m_V6_i_max = m_V5_i_max = i;
        }
        else if (pt.m_v >= V0() - (m_V7+m_V6+m_V5))
        {
            m_V5_i_max = i;
        }
        else
        {
            break;
        }
    }
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

void Process::exec()
{
    m_detail->first_h(m_v_parts);
    m_detail->next_h(1.0);
    //for (m_detail->first_h(m_v_parts); m_detail->isValid(); m_detail->next_h())
    {
        // ...
    }
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
