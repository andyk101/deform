#ifndef DETAIL_H
#define DETAIL_H
// -----------------------------------------------------------------------------------------------------------

#include <QtCore>
#include <andyk/core/core_fwd.h>

// -----------------------------------------------------------------------------------------------------------
// Material
// -----------------------------------------------------------------------------------------------------------
enum eDirection { eDeg0, eDeg45, eDegCount };
struct Material
{
    QString m_name;
    double m_R[eDegCount];
    double m_B;
    double m_m;
    double m_Omega;
    double m_U;
    double m_aa0;
    double m_aa1;
    double m_aa2;

public:
    static QSharedPointer<Material> find(QString name);

    Material(){}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
};

// -----------------------------------------------------------------------------------------------------------
// Detail - деталь при шаге h
// -----------------------------------------------------------------------------------------------------------
struct Geom;
class Detail
{
public:
    static QSharedPointer<Detail> find(QString name);

private:
    QString m_name;
    QSharedPointer<Material> m_material;
    QSharedPointer<Geom> m_geom[eDegCount];

public:
    // парсинг
    double m_r_0;
    double m_m_d;
    double m_r_kp;
    double m_r_km;
    double m_s_0;
    double m_z;
    double m_mu;

    // инициализация
    double m_r_2;
    double m_r_c;
    double m_r_max;
    double m_V0;
    double m_V1;

public:
    Detail() { init(); }
    Detail& operator= (const Detail& detail);
    Detail(const Detail& detail) { init(); *this = detail; }
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    QSharedPointer<Detail> clone() const;

    QString name() const { return m_name; }
    QSharedPointer<Material> material() const { return m_material; }
    void setMaterial(QSharedPointer<Material> material) { m_material = material; }
    QSharedPointer<Geom> geom(int direction) { return m_geom[direction]; }

    void init();
    void first(int v_parts);
    bool isValid() const;

    double h() const;
    void next_h(double dh);
    void calcContext(QSharedPointer<Geom> prevGeom[eDegCount], bool calc_s);
    void recalc_r_max();
};

// -----------------------------------------------------------------------------------------------------------
// Point - точка геометрии c шагом по объему v
// -----------------------------------------------------------------------------------------------------------
struct Geom;
struct GeomsPoint
{
    Geom* m_geom;

    double m_v;
    double m_r;
    double m_h;
    double m_alpha;
    double m_s;    

    double m_epsilon_phi;
    double m_epsilon_i;
    double m_sigma_s;
    double m_sigma_r;
    double m_sigma_phi;
    double m_s_expr;
    double m_omega_e;

    GeomsPoint() { init(0,0); }
    GeomsPoint(Geom* geom, double s) { init(geom,s); }
    void init(Geom* geom, double s);
};

// -----------------------------------------------------------------------------------------------------------
// Geom - геометрия детали в направлении d
// -----------------------------------------------------------------------------------------------------------
enum eV { eV7, eV6, eV5, eVCount };
struct Geom
{
    Material* m_material; // материал
    Detail* m_detail;     // деталь
    int m_direction;      // направление d
    QVector<GeomsPoint> m_points;
    GeomsPoint m_pt_rmax;

    double m_h, m_max_dh;

    double m_s_1;
    double m_V2;
    double m_V5;
    double m_V6;
    double m_V7;
    double m_V5_bound;
    double m_V6_bound;
    double m_V7_bound;
    double m_alpha;
    double m_AB;
    double m_r_1;
    double m_r_k;

    double m_V_eps;
    int m_V7_i_max;
    int m_V6_i_max;
    int m_V5_i_max;

    // fV6(alpha_xx)+fV5(alpha_xx)+fV2(alpha_xx)=v <-> a5*sin(x/2)^2*cos(x)+b5*x*cos(x)+a3*sin(x)^2+b3*sin(x)-v*cos(x)+c3=0
    //               fV5(alpha_xx)+fV2(alpha_xx)=v <-> a1*sin(x/2)^2*cos(x)+b1*x*cos(x)+a3*sin(x)^2+b3*sin(x)-v*cos(x)+c3=0
    //                               fV6(alpha_xx) <-> a4*sin(alpha_xx/2)^2+b4*alpha_xx;
    //                               fV5(alpha_xx) <-> (a3*sin(alpha_xx)^2+b3*sin(alpha_xx)+c3)/cos(alpha_xx);
    //                               fV2(alpha_xx) <-> a1*sin(alpha_xx/2)^2+b1*alpha_xx;
    double m_a1, m_b1;
    double m_a2, m_b2, m_c2, m_d2;
    double m_a3, m_b3, m_c3;
    double m_a4, m_b4;
    double m_a5, m_b5;

    bool m_valid;

    // материал
    double R(int direction) const { return m_material->m_R[direction]; }
    double R() const { return m_material->m_R[m_direction]; }
    double B() const { return m_material->m_B; }
    double m() const { return m_material->m_m; }
    double Omega() const { return m_material->m_Omega; }
    double U() const { return m_material->m_U; }
    double aa0() const { return m_material->m_aa0; }
    double aa1() const { return m_material->m_aa1; }
    double aa2() const { return m_material->m_aa2; }

    // деталь
    double r_0() const { return m_detail->m_r_0; }
    double m_d() const { return m_detail->m_m_d; }
    double r_kp() const { return m_detail->m_r_kp; }
    double r_km() const { return m_detail->m_r_km; }
    double s_0() const { return m_detail->m_s_0; }
    double z() const { return m_detail->m_z; }
    double mu() const { return m_detail->m_mu; }
    double r_2() const { return m_detail->m_r_2; }
    double r_c() const { return m_detail->m_r_c; }
    double r_max() const { return m_detail->m_r_max; }
    double V0() const { return m_detail->m_V0; }
    double V1() const { return m_detail->m_V1; }

    // дополнительные
    int eV_x(double v_xx) const;
    double AB_x(double alpha_xx) const;
    double h_x(double alpha_xx, double AB_xx) const;
    double V2_x(double alpha_xx) const;
    double V5_x(double alpha_xx, double AB_xx) const;
    double V6_x(double alpha_xx) const;
    double V7_x(double r_k_xx) const;

    bool isValid() const { return m_valid; }

    void calcPoint(double& r_x, double& h_x, double& alpha_x, double v_x) const;
    Geom(Detail* detail, int direction, int count, double dv);
    QSharedPointer<Geom> clone() const { return QSharedPointer<Geom>(new Geom(*this)); }
    void recalc_V_coeff();
    void recalc_max_dh();
    void recalc_r_max();

    void next_h(double dh);
    void calcContext(QSharedPointer<Geom> prevGeom[eDegCount], bool calc_s);
};

// -----------------------------------------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------------------------------------
class Process
{
public:
    static QSharedPointer<Process> find(QString name);
    struct GeomPair
    {
        QSharedPointer<Geom> m_geoms[eDegCount];

        GeomPair() {}
        GeomPair(QSharedPointer<Geom> geom0, QSharedPointer<Geom> geom45)
        {
            m_geoms[eDeg0] = geom0;
            m_geoms[eDeg45] = geom45;
        }
    };

private:
    QString m_name;
    QList< GeomPair > m_prev_geoms; // для отладки
    QSharedPointer<Detail> m_detail;
    GeomsPoint m_pt_r_max[eDegCount];
    int m_v_parts;
    bool m_calc_s;

public:
    Process() {}
    Process& operator= (const Process& proc);
    Process(const Process& proc) { *this = proc; }
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    QSharedPointer<Process> clone() const;

    QString name() const { return m_name; }
    QSharedPointer<Detail> detail() const { return m_detail; }
    void setDetail(QSharedPointer<Detail> detail) { m_detail = detail; }

    void exec();
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

// Возможно критерии тоже не надо клонировать???
// Для этого критерий не должен содержать Process.
// C клонированием гибче, т.к. возможно что-то будет особое для каждого критерия хрантиться???
class Criterion
{
public:
    static QSharedPointer<Criterion> find(QString name);

private:
    QString m_name;
    QSharedPointer<Process> m_process;
    CriterionType m_type;
    double m_method;

public:
    Criterion() {}
    Criterion& operator= (const Criterion& crit);
    Criterion(const Criterion& crit) { *this = crit; }
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    QSharedPointer<Criterion> clone() const;

    QString name() const { return m_name; }
    QSharedPointer<Process> process() const { return m_process; }
    void setProcess(QSharedPointer<Process> proc) { m_process = proc; }
};

// -----------------------------------------------------------------------------------------------------------
// Curve
// -----------------------------------------------------------------------------------------------------------
struct CurveArgTypes : public Enum
{
    enum { e_r_km, e_r_kp, e_mu };

    CurveArgTypes() : Enum("CurveArgTypes")
    {
        *this << "r_km" << "r_kp" << "mu";
    }
};
typedef EnumItem<CurveArgTypes> CurveArgType;

class Curve
{
public:
    static QSharedPointer<Curve> find(QString name);

private:
    QString m_name;
    CurveArgType m_argType;
    double m_argStart;
    double m_argEnd;
    double m_argStep;

public:
    Curve() {}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
};

// -----------------------------------------------------------------------------------------------------------
// Plots
// -----------------------------------------------------------------------------------------------------------
struct PlotsCollections : public Enum
{
    enum { eDetail, eProcess, eCurve, eMaterial, eCriterion, eCount };

    PlotsCollections() : Enum("PlotsCollections")
    {
        *this << "Detail" << "Process" << "Curve" << "Material" << "Criterion";
    }
};
typedef EnumItem<PlotsCollections> PlotsCollection;

class Plots
{
public:
    static QSharedPointer<Plots> find(QString name);
    static QString cartesianName(QString strTemplate, QStringList strCartesian);

private:
    QString m_name;
    QString m_plotName;
    QStringList m_collections[PlotsCollections::eCount];

public:
    Plots() {}
    void addString(QDomElement& elem, QStringList* lst);
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    void create();
};

// -----------------------------------------------------------------------------------------------------------
// Plot
// -----------------------------------------------------------------------------------------------------------
struct Plot
{
public:
    static QSharedPointer<Plot> find(QString name);

private:
    QString m_name;
    QSharedPointer<Criterion> m_criterion;
    QSharedPointer<Curve> m_curve;
    bool m_empty;

public:
    Plot() {}

    bool isEmpty() const { return m_empty; }
    void reset(QString name, QStringList strCartesian);

    QString name() const { return m_name; }
    QSharedPointer<Criterion> criterion() const { return m_criterion; }
    QSharedPointer<Curve> curve() const { return m_curve; }

    void calculate();
};

// -----------------------------------------------------------------------------------------------------------
// Layouts
// -----------------------------------------------------------------------------------------------------------
class Layouts
{
public:
    static QSharedPointer<Layouts> find(QString name);

private:
    QString m_name;
    QString m_layoutName;
    QStringList m_plots;
    QStringList m_collections[PlotsCollections::eCount];

public:
    Layouts() {}
    void addString(QDomElement& elem, QStringList* lst);
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    void create();
};

// -----------------------------------------------------------------------------------------------------------
// Layout
// -----------------------------------------------------------------------------------------------------------
struct Layout
{
public:
    static QSharedPointer<Layout> find(QString name);

private:
    QString m_name;
    QList<QSharedPointer<Plot> > m_plots;
    bool m_empty;

public:
    Layout() {}

    bool isEmpty() const { return m_empty; }
    void reset(QString name, QStringList plotNames);

    QString name() const { return m_name; }
    const QList<QSharedPointer<Plot> >& plots() const { return m_plots; }

    void create();
};


// -----------------------------------------------------------------------------------------------------------
// Point_h_v - точка в направлении d, при шаге h, и объеме v
// -----------------------------------------------------------------------------------------------------------
//double d() const;
//double h() const;

// -----------------------------------------------------------------------------------------------------------
#endif // DETAIL_H
