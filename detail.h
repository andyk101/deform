#ifndef DETAIL_H
#define DETAIL_H
// -----------------------------------------------------------------------------------------------------------

#include <andyk/core/core_fwd.h>
#include <andyk/serialization/serialization.h>

// -----------------------------------------------------------------------------------------------------------
// Material
// -----------------------------------------------------------------------------------------------------------
using namespace boost::serialization;
struct Material
{
    std::string m_name;
    double m_R0;
    double m_R45;
    double m_B;
    double m_m;
    double m_Omega;
    double m_U;
    double m_aa0;
    double m_aa1;
    double m_aa2;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
        ar & make_attr_nvp("Name", m_name);
        ar & make_attr_nvp("R0", m_R0);
        ar & make_attr_nvp("R45", m_R45);
        ar & make_attr_nvp("B", m_B);
        ar & make_attr_nvp("m", m_m);
        ar & make_attr_nvp("Omega", m_Omega);
        ar & make_attr_nvp("U", m_U);
        ar & make_attr_nvp("aa0", m_aa0);
        ar & make_attr_nvp("aa1", m_aa1);
        ar & make_attr_nvp("aa2", m_aa2);
    }

public:
    static boost::shared_ptr<Material> find(std::string name);

    Material(){}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
};

/*
// -----------------------------------------------------------------------------------------------------------
// Point_h_v - точка в направлении d, при шаге h, и объеме v
// -----------------------------------------------------------------------------------------------------------
struct Geom;
struct Point;
{
    Geom* m_geom;

    double m_v;
    double m_r;
    double m_h;
    double m_s;

    double m_epsilon_phi;
    double m_epsilon_i;
    double m_sigma_s;
    double m_sigma_r;
    double m_sigma_phi;
    double m_s_expr;
    double m_omega_e;

    // ...

    Point_d_h_v()
    {
    }

    double d() const;
    double h() const;
};

// -----------------------------------------------------------------------------------------------------------
// Geom - геометрия в направлении d, при шаге h
// -----------------------------------------------------------------------------------------------------------
struct Geom
{
    Geom* m_pGeom_prev_h;
    QVector<> m_points;
    double m_h;

    double m_s_1_0;
    double m_s_1_45;
    double m_h;
    double m_V2;
    double m_V5;
    double m_V6;
    double m_V7;
    double m_alpha;
    double m_AB;
    double m_r_1;
    double m_r_k;


public:
    Geom();


};

// -----------------------------------------------------------------------------------------------------------
// Direction - направление
// -----------------------------------------------------------------------------------------------------------
struct Direction
{
public:

private:
    Geom m_geom[eDegSize];

public:
    Direction()
    {
    }
};
*/

// -----------------------------------------------------------------------------------------------------------
// Detail
// -----------------------------------------------------------------------------------------------------------
class Detail
{
public:
    static boost::shared_ptr<Detail> find(std::string name);

private:
    friend class boost::serialization::access;

    std::string m_name;
    boost::shared_ptr<Material> m_material;

    double m_r_0;
    double m_m_d;
    double m_r_kp;
    double m_r_km;
    double m_s_0;
    double m_z;
    double m_mu;

    double m_r_2;
    double m_r_c;
    double m_r_max;
    double m_V0;
    double m_V1;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
        ar & make_attr_nvp("Name", m_name);
        ar & make_attr_nvp("r_0", m_r_0);
        ar & make_attr_nvp("m_d", m_m_d);
        ar & make_attr_nvp("r_kp", m_r_kp);
        ar & make_attr_nvp("r_km", m_r_km);
        ar & make_attr_nvp("s_0", m_s_0);
        ar & make_attr_nvp("z", m_z);
        ar & make_attr_nvp("mu", m_mu);
        ar & make_nvp("material", m_material);
    }

    //QVector<Direction> m_direction;
public:
    Detail() {}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    boost::shared_ptr<Detail> clone() const;
};

// -----------------------------------------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------------------------------------
class Process
{
public:
    static boost::shared_ptr<Process> find(std::string name);

private:
    friend class boost::serialization::access;

    std::string m_name;
    boost::shared_ptr<Detail> m_detail;

    int m_v_parts;
    bool m_calc_s;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
        ar & make_attr_nvp("Name", m_name);
        ar & make_attr_nvp("v_parts", m_v_parts);
        ar & make_attr_nvp("calc_s", m_calc_s);
        ar & make_nvp("detail", m_detail);
    }

public:
    Process() {}
    Process& operator= (const Process& proc);
    Process(const Process& proc) { *this = proc; }
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    boost::shared_ptr<Process> clone() const;

    //void step();
};

// -----------------------------------------------------------------------------------------------------------
// Critarion
// -----------------------------------------------------------------------------------------------------------
struct CriterionTypes : public Enum
{
    enum { eFenom, eSigmaRMax, eLocal };

    CriterionTypes() : Enum("CriterionTypes")
    {
        *this << "Fenom" << "SigmaRMax" << "Local";
    }
};
typedef EnumItem<CriterionTypes> CriterionType;

class Criterion
{
public:
    static boost::shared_ptr<Criterion> find(std::string name);

private:
    friend class boost::serialization::access;

    std::string m_name;
    boost::shared_ptr<Process> m_process;
    CriterionType m_type;
    double m_method;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
        ar & make_attr_nvp("Name", m_name);

        // wrapper_traits!

//        std::string strType;
//        if (Archive::is_saving::value)
//        {
//            strType = sm_criterionType.item(m_valType).toStdString();
//            ar << make_attr_nvp("Type", strType);
//        }
//        else
//        {
//            ar >> make_attr_nvp("Type", strType);
//            m_valType = sm_criterionType.item(QString::fromStdString(strType));
//        }

        ar & make_attr_nvp("Method", m_method);
        ar & make_nvp("process", m_process);
    }

public:
    Criterion() {}
    Criterion& operator= (const Criterion& crit);
    Criterion(const Criterion& crit) { *this = crit; }
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    boost::shared_ptr<Criterion> clone() const;
};

// -----------------------------------------------------------------------------------------------------------
// Curve - m_d от аргумента
// -----------------------------------------------------------------------------------------------------------
struct Curve
{
};


// -----------------------------------------------------------------------------------------------------------
// Plot - график
// -----------------------------------------------------------------------------------------------------------
struct Plot
{
    QList<QWeakPointer<Detail> > m_detal;
    QList<QWeakPointer<Curve> > m_curve;
    QList<QWeakPointer<Material> > m_material;
    QList<QWeakPointer<Criterion> > m_critarion;

    Plot(QWeakPointer<Detail> detal, QWeakPointer<Detail> curve, QWeakPointer<Detail> material,
         QWeakPointer<Criterion> critarion);
};

// -----------------------------------------------------------------------------------------------------------
// Point_h_v - точка в направлении d, при шаге h, и объеме v
// -----------------------------------------------------------------------------------------------------------
//double d() const;
//double h() const;

// -----------------------------------------------------------------------------------------------------------
#endif // DETAIL_H
