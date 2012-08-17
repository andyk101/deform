#include "detail.h"
#include "deform_app.h"
#include <boost/math/special_functions/sign.hpp>

// ----------------------------------------------------------------------------------------------------------
// Решение квадратного уравнения a*x^2+b*x+c=0, если a,b,c - числа с плав. точкой.
// ----------------------------------------------------------------------------------------------------------
int SquareSolve(double& x1, double& x2, double a, double b, double c)
{
    if (a == 0)
    {
        // c=0
        if (b == 0)
            return false;

        // b*x+c=0
        x1 = -c/b;
        return 1;
    }

    // x^2+b*x+c=0
    b /= a;
    c /= a;

    double d = b*b-4*c;
    if (d < 0)
        return 0;

    d = sqrt(d)/2;
    b = -b/2;

    if (d == 0)
    {
        x1 = b;
        return 1;
    }

    x1 = b-d;
    x2 = b+d;
    return 2;
}

// ----------------------------------------------------------------------------------------------------------
// Решение кубич. уравнения a*x^3+b*x^2+c*x+d=0, если a,b,c,d - числа с плав. точкой.
// ----------------------------------------------------------------------------------------------------------
int CubeSolve(double& x1, double& x2, double& x3, double a, double b, double c, double d)
{
    if (a == 0)
    {
        return SquareSolve(x1, x2, b, c, d);
    }
    double p = -pow(b,2)/(3*pow(a,2))+c/a;
    double q = 2*pow(b,3)/(27*pow(a,3))-b*c/(3*pow(a,2))+d/a;
    double Q = pow(p/3,3)+pow(q/2,2);
    double R = boost::math::sign(q)*sqrt(fabs(p)/3);

    // 1 корень
    double phi = 0;
    if (p > 0)
    {
        phi = asinh(q/(2*pow(R,3)));
        x1 = -2*R*sinh(phi/3)-b/(3*a);
        return 1;
    }
    if (Q > 0)
    {
        phi = acosh(q/(2*pow(R,3)));
        x1 = -2*R*cosh(phi/3)-b/(3*a);
        return 1;
    }

    // 3 корня
    phi = acos(q/(2*pow(R,3)));
    x1 = -2*R*cos(phi/3)-b/(3*a);
    x2 = -2*R*cos(phi/3+2*M_PI/3)-b/(3*a);
    x3 = -2*R*cos(phi/3+4*M_PI/3)-b/(3*a);
    return 3;
}

// -----------------------------------------------------------------------------------------------------------
// Нахождение корня func.f(x)=0 методом Итераций/Ньютона
// ----------------------------------------------------------------------------------------------------------
template<class Func> bool NewtonSolve(Func func, double& x, double x0, double eps)
{
    double x1 = func.phi(x0);
    double x2 = func.phi(x1);
    double diff2 = fabs(x2 - x1);
    double diff1 = 2*diff2;

    while ( diff2 < diff1 && diff2 > 0 )
    {
        x1 = x2;
        diff1 = diff2;
        x2 = func.phi(x1);
        diff2 = fabs(x2 - x1);
    }
    if (diff2 > eps)
        return false;

    x = x2;
    return true;
}

// -----------------------------------------------------------------------------------------------------------
// Сужение интервала поиска корня методом Дихотомии и последующее нахождение корня методом Ньютона
// ----------------------------------------------------------------------------------------------------------
template<class Func> bool DichotomNewtonSolve(Func func, double& x, double x1, double x2, double eps)
{
    if (x1 > x2)
        throw ErrorBase::create(QString("Expected x1<=x2 (x1=%1, x2=%2)")
            .arg(x1).arg(x2));

    double f1 = func.f(x1);
    double f2 = func.f(x2);

    // test Dichotom condition
    if (boost::math::sign(f1*f2) > 0)
        return false;

    int splits_count = 0;
    while ( fabs(x2-x1) > eps )
    {
        // next x
        double x3 = (x1+x2)/2;

        // test Newton condition
        if (splits_count >= 4)
        {
            if (fabs(func.d2f(x1)) < pow(func.df(x1),2)/2 &&
                fabs(func.d2f(x2)) < pow(func.df(x2),2)/2)
            {
                // Newton method
                return NewtonSolve(func, x, x3, eps);
            }
            splits_count = 0;
        }

        // split interval
        double f3 = func.f(x3);
        if (boost::math::sign(f1*f3) <= 0)
        {
            x2 = x3;
            f2 = f3;
        }
        else
        {
            x1 = x3;
            f1 = f3;
        }
        splits_count++;
    }

    // Dichotom method
    return (x1+x2)/2;
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

void Detail::first(int v_parts)
{
    m_r_2 = m_m_d*m_r_0 - m_z/2;
    m_r_c = m_r_km + m_z + m_r_2;
    m_r_max = m_m_d*m_r_0;
    m_V0 = M_PI*m_s_0*pow(m_r_0, 2);
    m_V1 = M_PI*m_s_0*pow(m_r_2-m_r_kp, 2);

//    double V7 = M_PI*m_s_0*( pow(m_r_0,2) - pow(m_r_c,2) );
//    int V7_parts = (int)round(V7 / (m_V0-m_V1) * v_parts);
//    double dv = V7 / V7_parts;
//    int count = (int)floor((m_V0-m_V1) / dv) + 1;
    double dv = m_V0 / v_parts;
    int count = (int)floor((m_V0-m_V1) / dv) + 1;

    m_geom[eDeg0] = QSharedPointer<Geom>(new Geom(this, eDeg0, count, dv));
    m_geom[eDeg45] = QSharedPointer<Geom>(new Geom(this, eDeg45, count, dv));
}

bool Detail::isValid() const
{
    return m_geom[eDeg0]->isValid() && m_geom[eDeg45]->isValid();
}

void Detail::next(double dh)
{
    dh = qMin(qMin(dh, m_geom[eDeg0]->m_max_dh), m_geom[eDeg45]->m_max_dh);
    m_geom[eDeg0]->next(dh);
    //m_geom[eDeg45]->next(dh);
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
    m_v = m_r = m_h = m_alpha = m_epsilon_phi = m_epsilon_i =
          m_sigma_s = m_sigma_r = m_sigma_phi = m_s_expr = m_omega_e = 0;
}

// -----------------------------------------------------------------------------------------------------------
// Решение уравнения a4*sin(x/2)^2+b4*x-v=0
// ----------------------------------------------------------------------------------------------------------
// fV6(alpha_xx) <-> a4*sin(alpha_xx/2)^2+b4*alpha_xx
// ----------------------------------------------------------------------------------------------------------
struct FuncV6
{
    double a4, b4, v;

    double f(double x) const
    {
        return a4*pow(sin(x/2),2)+b4*x-v;
    }
    double df(double x) const
    {
        return a4/2*sin(x)+b4;
    }
    double phi(double x) const
    {
        return x - f(x)/df(x);
    }
};

// -----------------------------------------------------------------------------------------------------------
// Решение уравнения ax*sin(x/2)^2*cos(x)+bx*x*cos(x)+a3*sin(x)^2+b3*sin(x)-v*cos(x)+c3=0
// -----------------------------------------------------------------------------------------------------------
// fV6(alpha_xx)+fV5(alpha_xx)+fV2(alpha_xx)=v <-> a5*sin(x/2)^2*cos(x)+b5*x*cos(x)+a3*sin(x)^2+b3*sin(x)-v*cos(x)+c3=0
//               fV5(alpha_xx)+fV2(alpha_xx)=v <-> a1*sin(x/2)^2*cos(x)+b1*x*cos(x)+a3*sin(x)^2+b3*sin(x)-v*cos(x)+c3=0
// ----------------------------------------------------------------------------------------------------------
struct FuncSumV
{
    double ax, bx, a3, b3, c3, v;

    FuncSumV(){}
    FuncSumV(double ax_, double bx_, double a3_, double b3_, double c3_, double v_)
        : ax(ax_), bx(bx_), a3(a3_), b3(b3_), c3(c3_), v(v_)
    {
    }
    double f(double x) const
    {
        return ax*pow(sin(x/2),2)*cos(x)+bx*x*cos(x)+a3*pow(sin(x),2)+b3*sin(x)-v*cos(x)+c3;
    }
    double df(double x) const
    {
        return -ax*pow(sin(x/2),2)*sin(x)-bx*x*sin(x)+(ax/4+a3)*sin(2*x)+(bx+b3)*cos(x)+v*sin(x);
    }
    double d2f(double x) const
    {
        return -ax*pow(sin(x/2),2)*cos(x)-bx*x*cos(x)+(-3*ax/2-4*a3)*pow(sin(x),2)-(2*bx+b3)*sin(x)+v*cos(x)+ax/2+2*a3;
    }
    double phi(double x) const
    {
        return x - f(x)/df(x);
    }
};

// -----------------------------------------------------------------------------------------------------------
// Geom
// -----------------------------------------------------------------------------------------------------------
Geom::Geom(Detail* detail, int direction, int count, double dv)
    : m_material(&*detail->material()), m_detail(detail), m_direction(direction),
      m_points(count, GeomsPoint(this, s_0())),
      m_valid(false)
{
    m_h = m_max_dh = 0;
    m_s_1 = s_0();
    m_V2 = m_V6 = 0;
    m_V5 = M_PI*s_0()*( pow(r_c(),2) - pow(r_2()-r_kp(),2) );
    m_V7 = M_PI*s_0()*( pow(r_0(),2) - pow(r_c(),2) );
    m_V5_bound = V0()-m_V7-m_V6-m_V5;
    m_V6_bound = V0()-m_V7-m_V6;
    m_V7_bound = V0()-m_V7;
    m_alpha = 0;
    m_AB = r_c() - (r_2() - r_kp());
    m_r_1 = r_c();
    m_r_k = r_0();

    // Т.к. объемы считаются не точно, точка переходит на следующий участок чуть раньше.
    m_V_eps = 0.00001;
    m_V7_i_max = m_V6_i_max = m_V5_i_max = -1;
    for (int i = 0; i < m_points.count(); i++)
    {
        GeomsPoint& pt = m_points[i];
        double v = V0() - i*dv;
        switch (eV_x(v))
        {
            case eV7:
                m_V7_i_max = m_V6_i_max = m_V5_i_max = i;
                break;
            case eV5:
                m_V5_i_max = i;
                break;
            default:
                goto Break;
        }

        pt.m_v = v;
        pt.m_r = (i == 0) ? r_0() :
                            sqrt(pow(m_points[i-1].m_r,2) - dv/(M_PI*s_0()));
    }
Break:

    recalc_V_coeff();
    recalc_max_dh();
    m_valid = true;
}

void Geom::recalc_V_coeff()
{
    m_a1 = 4*M_PI/3*s_0()*(3*pow(r_kp(),2)+3*r_kp()*s_0()+pow(s_0(),2));
    m_b1 = M_PI*s_0()*(2*r_kp()+s_0())*(r_2()-r_kp());

    m_a2 = r_km()+r_kp()+s_0();
    m_b2 = r_km()*s_0()+2*r_km()*m_s_1-2*r_kp()*s_0()-r_kp()*m_s_1-pow(s_0(),2)+pow(m_s_1,2);
    m_c2 = -r_c()+r_2()-r_kp();
    m_d2 = -r_c()*s_0()-2*r_c()*m_s_1-2*r_2()*s_0()-r_2()*m_s_1+2*r_kp()*s_0()+m_s_1*r_kp();

    m_a3 = M_PI/3*m_a2*m_b2;
    m_b3 = M_PI/3*(m_a2*m_d2+m_c2*m_b2);
    m_c3 = M_PI/3*m_c2*m_d2;

    m_a4 = -4*M_PI/3*m_s_1*(3*pow(r_km(),2)+3*r_km()*m_s_1+pow(m_s_1,2));
    m_b4 = M_PI*m_s_1*(2*r_km()+m_s_1)*r_c();

    m_a5 = m_a1+m_a4;
    m_b5 = m_b1+m_b4;
}

void Geom::recalc_max_dh()
{
    if (m_V6_i_max < 0)
    {
        m_max_dh = 0;
        return;
    }

    m_max_dh = std::numeric_limits<double>::max();
    if (m_V7_i_max >= 0)
    {
        double dv = V7_x(m_points[m_V7_i_max].m_r);

        double alpha_xx = 0;
        FuncSumV func(m_a5, m_b5, m_a3, m_b3, m_c3, m_V2+m_V5+m_V6+dv);
        if(!DichotomNewtonSolve(func, alpha_xx, m_alpha, M_PI/2, 0.000000000001))
            throw ErrorBase::create(QString("Cannot find alpha_x for dv=%1").arg(dv));

        double AB = AB_x(alpha_xx);
        double vx = V6_x(alpha_xx)+V5_x(alpha_xx, AB)+V2_x(alpha_xx);
        double eps = fabs(func.v - vx);
        double max_dh = h_x(alpha_xx, AB) - m_h;
        m_max_dh = qMin(m_max_dh, max_dh);
    }

    if (m_V6_i_max > m_V7_i_max)
    {
        double dv = m_V6 - V6_x(m_points[m_V6_i_max].m_alpha);

        double alpha_xx = 0;
        FuncSumV func(m_a1, m_b1, m_a3, m_b3, m_c3, m_V2+m_V5+dv);
        if(!DichotomNewtonSolve(func, alpha_xx, m_alpha, M_PI/2, 0.000000000001))
            throw ErrorBase::create(QString("Cannot find alpha_x for dv=%1").arg(dv));

        double AB = AB_x(alpha_xx);
        double vx = V5_x(alpha_xx, AB)+V2_x(alpha_xx);
        double eps = fabs(func.v - vx);
        double max_dh = h_x(alpha_xx, AB) - m_h;
        m_max_dh = qMin(m_max_dh, max_dh);
    }
}

// -----------------------------------------------------------------------------------------------------------
int Geom::eV_x(double v_xx) const
{
    return           (v_xx > V0()) ? eVCount :
     (v_xx > m_V7_bound + m_V_eps) ? eV7 :
     (v_xx > m_V6_bound + m_V_eps) ? eV6 :
     (v_xx > m_V5_bound + m_V_eps) ? eV5 :
                                     eVCount;
}

double Geom::AB_x(double alpha_xx) const
{
    return r_c()-r_2()+r_kp()-(sin(alpha_xx)*(r_km()+r_kp()+s_0()));
}

double Geom::h_x(double alpha_xx, double AB_xx) const
{
    return AB_xx*tan(alpha_xx)+(1-cos(alpha_xx))*(r_km()+r_kp()+s_0());
}

double Geom::V2_x(double alpha_xx) const
{
    return 4*M_PI/3*s_0()*pow(sin(alpha_xx/2),2)*(3*pow(r_kp(),2) + 3*r_kp()*s_0() + pow(s_0(),2)) +
        M_PI*alpha_xx*s_0()*(2*r_kp()+s_0())*(r_2()-r_kp());
}

double Geom::V5_x(double alpha_xx, double AB_xx) const
{
    return M_PI*AB_xx/(3*cos(alpha_xx)) * ( AB_xx*(s_0()+2*m_s_1) + 3*(s_0()+m_s_1)*(r_2()-r_kp()) +
        sin(alpha_xx)*((s_0()+m_s_1)*(3*r_kp()+2*s_0())-pow(m_s_1,2)) );
}

double Geom::V6_x(double alpha_xx) const
{
    return M_PI*alpha_xx*m_s_1*(2*r_km()+m_s_1)*r_c() - 4.0/3*M_PI*m_s_1*pow(sin(alpha_xx/2),2)*
        (3*pow(r_km(),2)+3*r_km()*m_s_1+pow(m_s_1,2));
}

double Geom::V7_x(double r_k_xx) const
{
    return M_PI*m_s_1*(pow(r_k_xx,2) - pow(r_c(),2));
}

// -----------------------------------------------------------------------------------------------------------
void Geom::calcPoint(double& r_x, double& h_z, double& alpha_x, double v_x) const
{
    if (m_V7_bound < v_x && v_x <= m_V7_bound + m_V_eps)
        v_x = m_V7_bound;
    else if (m_V6_bound < v_x && v_x <= m_V6_bound + m_V_eps)
        v_x = m_V6_bound;
    else if (m_V5_bound < v_x && v_x <= m_V5_bound + m_V_eps)
        v_x = m_V5_bound;

    switch(eV_x(v_x))
    {
        case eVCount:
            throw ErrorBase::create(QString("v_x must be in range [%1, %2], but actual is %3")
                .arg(m_V5_bound).arg(V0()).arg(v_x));

        case eV7:
            alpha_x = 0;
            r_x = sqrt((v_x-V0()+m_V7)/(M_PI*m_s_1)+pow(r_c(),2));
            h_z = 0;
            break;

        case eV6:
        {
            FuncV6 func;
            func.a4 = m_a4;
            func.b4 = m_b4;
            func.v = V0()-m_V7-v_x;

            if (!NewtonSolve(func, alpha_x, M_PI/4, 0.000000000001))
                throw ErrorBase::create(QString("Cannot find alpha_x for v_x=%1").arg(v_x));

            r_x = r_c()-(r_km()+s_0()/2)*sin(alpha_x);
            h_z = (r_km()+s_0()/2)*(1-cos(alpha_x));
            break;
        }
        case eV5:
        {
            // fV5(AB_xx,fs_1(AB_xx))=v-V1-V2) <-> { a1*AB_xx^3+b1*AB_xx^2+c1*AB_xx+d1=0, 0<=AB_xx<=AB }
            double a = M_PI*(m_s_1-s_0())*(-2*m_AB+sin(m_alpha)*(m_s_1-s_0()));
            double b = -3*M_PI*m_AB*(m_AB*s_0()+(m_s_1-s_0())*((r_2()-r_kp())+sin(m_alpha)*r_kp()));
            double c = 3*M_PI*pow(m_AB,2)*s_0()*(-2*r_2()+2*r_kp()-sin(m_alpha)*(2*r_kp()+s_0()));
            double d = 3*cos(m_alpha)*pow(m_AB,2)*(v_x-V1()-m_V2);

            double AB_x, x1,x2,x3;
            int roots = CubeSolve(x1,x2,x3, a,b,c,d);
            if (roots>=1 && 0<=x1 && x1<=m_AB)
                AB_x = x1;
            else if (roots>=2 && 0<=x2 && x2<=m_AB)
                AB_x = x2;
            else if (roots>=3 && 0<=x3 && x3<=m_AB)
                AB_x = x3;
            else
                throw ErrorBase::create(QString("Cannot find AB_x for v_x=%1").arg(v_x));

            alpha_x = 0;
            r_x = m_r_1 - m_AB + AB_x;
            h_z = m_h - (r_kp()+s_0()/2)*(1-cos(m_alpha)) - tan(m_alpha)*AB_x;
            break;
        }
    }
}

void Geom::next(double dh)
{
    if (!isValid())
        return;

    // m_h
    dh = qMin(dh, m_max_dh);
    if (dh <= 0)
    {
        m_valid = false;
        return;
    }
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

        double x, x1, x2;
        int roots = SquareSolve(x1,x2, a2,b2,c2);
        if (roots>=1 && x1>=0 && (c1-x1*a1)/b1>=0)
            x = x1;
        else if (roots>=2 && x2>=0 && (c1-x2*a1)/b1>=0)
            x = x2;
        else
        {
            m_valid = false;
            return;
        }
        m_alpha = asin(x);
    }
    //m_alpha = M_PI / 4;
    //m_s_1 = s_0()*1.3;
    m_AB = r_c() - r_2() + r_kp() - sin(m_alpha)*(r_km()+r_kp()+s_0());
    m_r_1 = r_c() - (r_km()+s_0()/2)*sin(m_alpha);

    // -------------------------------------------------------------------------------------------------------
    // (s_0, s_1, alpha, AB, r_1) -> (V2, V5, V6, V7, V5_bound, V6_bound, V7_bound)
    m_V2 = V2_x(m_alpha);
    m_V5 = V5_x(m_alpha, m_AB);
    m_V6 = V6_x(m_alpha);
    if (V0() < V1()+m_V2+m_V5)
    {
        m_V5 = V0() - (V1()+m_V2);
        m_V6 = 0;
        m_V7 = 0;
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
    m_V5_bound = V0()-m_V7-m_V6-m_V5;
    m_V6_bound = V0()-m_V7-m_V6;
    m_V7_bound = V0()-m_V7;

    // -------------------------------------------------------------------------------------------------------
    // (points_prev) -> (r_k, points, V7_i_max, V6_i_max, V5_i_max)
    int V5_i_max_prev = m_V5_i_max;
    m_V7_i_max = m_V6_i_max = m_V5_i_max = -1;
    for (int i = 0; i <= V5_i_max_prev; i++)
    {
        GeomsPoint& pt = m_points[i];
        switch (eV_x(pt.m_v))
        {
            case eV7:
                m_V7_i_max = m_V6_i_max = m_V5_i_max = i;
                break;
            case eV6:
                m_V6_i_max = m_V5_i_max = i;
                break;
            case eV5:
                m_V5_i_max = i;
                break;
            case eVCount:
                goto Break;
        }

        calcPoint(pt.m_r, pt.m_h, pt.m_alpha, pt.m_v);

        if (i == 0)
            m_r_k = pt.m_r;
    }
Break:

    recalc_V_coeff();
    recalc_max_dh();
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
    for (m_detail->first(m_v_parts); m_detail->isValid(); m_detail->next(27))
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
