#ifndef DETAIL_H
#define DETAIL_H
// -----------------------------------------------------------------------------------------------------------

#include <QtCore>
#include <andyk/core/core_fwd.h>

// -----------------------------------------------------------------------------------------------------------
// Material
// -----------------------------------------------------------------------------------------------------------
struct Material
{
    QString m_name;
    double m_R0;
    double m_R45;
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
enum eDirection { eDeg0, eDeg45, eDegCount };
struct Geom;
class Detail
{
public:
    static QSharedPointer<Detail> find(QString name);

private:
    QString m_name;
    QSharedPointer<Material> m_material;
    QVector< QSharedPointer<Geom> > m_geom;

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
    int m_v_parts;
    double m_r_2;
    double m_r_c;
    double m_r_max;
    double m_V0;
    double m_V1;
    double m_dv;

public:
    Detail() : m_geom(eDegCount) {}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    QSharedPointer<Detail> clone() const;

    QString name() const { return m_name; }
    QSharedPointer<Material> material() const { return m_material; }
    void setMaterial(QSharedPointer<Material> material) { m_material = material; }
    QSharedPointer<Geom> geom(int direction) { return m_geom[direction]; }

    void first_h(int v_parts);
    bool isValid() const;
    void next_h(double dh);
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
    double m_s;

    double m_epsilon_phi;
    double m_epsilon_i;
    double m_sigma_s;
    double m_sigma_r;
    double m_sigma_phi;
    double m_s_expr;
    double m_omega_e;

    GeomsPoint();
    GeomsPoint(Geom* geom, double s);
};

// -----------------------------------------------------------------------------------------------------------
// Geom - геометрия детали в направлении d
// -----------------------------------------------------------------------------------------------------------
enum eBounds { eBoundV7, eBoundV6, eBoundRmax, eBoundV5, eBoundCount };
struct Geom
{
    Material* m_material; // материал
    Detail* m_detail;     // деталь
    int m_direction;      // направление d
    QVector<GeomsPoint> m_points;
    QVector<GeomsPoint> m_bounds;

    double m_h;
    double m_s_1;
    double m_V2;
    double m_V5;
    double m_V6;
    double m_V7;
    double m_alpha;
    double m_AB;
    double m_r_1;
    double m_r_k;

    int m_V7_i_max;
    int m_V6_i_max;
    int m_V5_i_max;

    bool m_valid;

    // материал
    double R0() const { return m_material->m_R0; }
    double R45() const { return m_material->m_R45; }
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
    int v_parts() const { return m_detail->m_v_parts; }
    double r_2() const { return m_detail->m_r_2; }
    double r_c() const { return m_detail->m_r_c; }
    double r_max() const { return m_detail->m_r_max; }
    double V0() const { return m_detail->m_V0; }
    double V1() const { return m_detail->m_V1; }
    double dv() const { return m_detail->m_dv; }

    // дополнительные
    double s_1_x(double AB_x) const;
    double V2_x(double alpha_x) const;
    double V5_x(double AB_x) const;
    double V6_x(double alpha_x) const;
    double V7_x(double r_k_x) const;
    bool isValid() const { return m_valid; }

    void calcPoint(double& r_x, double& h_x, double v_x) const;
    Geom(Detail* detail, int direction);
    void next_h(double dh);
};

// -----------------------------------------------------------------------------------------------------------
// Process
// -----------------------------------------------------------------------------------------------------------
class Process
{
public:
    static QSharedPointer<Process> find(QString name);

private:
    QString m_name;
    QSharedPointer<Detail> m_detail;

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
