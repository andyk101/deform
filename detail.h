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
    static QSharedPointer<Detail> find(QString name);

private:
    QString m_name;
    QSharedPointer<Material> m_material;

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

    //QVector<Direction> m_direction;
public:
    Detail() {}
    void parse(XmlParser& parser, QDomElement& element, bool optional = false);
    QSharedPointer<Detail> clone() const;

    QString name() const { return m_name; }
    QSharedPointer<Material> material() const { return m_material; }
    void setMaterial(QSharedPointer<Material> material) { m_material = material; }
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

    void start() {}
    void step() {}
    void exec() {}
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
