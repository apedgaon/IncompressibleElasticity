#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <Eigen/Dense>

class LinearLagrange
{
public:
    // 1D Element
    static double N(unsigned int index, double x);
    static double dNdx(unsigned int index, double x);

    // 2D Quadrilateral
    static double N_quad(unsigned int index, double xi, double eta);
    static double dNdxi_quad(unsigned int index, unsigned int didx, double xi, double eta);

    // 2D Triangle
    static double N_tri(unsigned int index, double xi, double eta);
    static double dNdxi_tri(unsigned int index, double xi, double eta);
    static double dNdeta_tri(unsigned int index, double xi, double eta);

    // 3D Hexahedral
    static double N_hex(unsigned int index, double xi, double eta, double kappa);
    static double dNdxi_hex(unsigned int index, unsigned int diff_idx, double xi, double eta, double kappa);
};

class QuadraticLagrange
{
public:
    // 1D Element
    static double N(unsigned int index, double x);
    static double dNdx(unsigned int index, double x);

    // 2D Quadrilateral
    static double N_quad(unsigned int index, double xi, double eta);
    static double dNdxi_quad(unsigned int index, unsigned int diff_idx, double xi, double eta);
};

class GsQuad
{
public:
    static double quadrature_pt(unsigned int num_pts, unsigned int index);
    static double weights(unsigned int num_pts, unsigned int index);
};

namespace Eigen {
    typedef Eigen::Matrix<unsigned int, -1, 1> VectorXui;
}
typedef std::vector<double> coords;
typedef std::vector<unsigned int> connectivity;

enum class MeshType
{
    linear = 0,
    quadrilateral = 1,
    triangular = 2,
    hexahedral = 3,
    tetrahedral = 4
};

enum class ShapeFcnType
{
    constant = 0,
    lagrange_linear = 1,
    lagrange_quadratic = 2,
    lagrange_cubic = 3
};

enum class QuadratureType
{
    one_point = 1,
    two_point = 2,
    three_point = 3
};

enum class AnalysisType
{
    _static = 1,
    _transient = 2,
    _modal = 3,
    _eigenbuckling = 4
};

template<unsigned int dim>
struct Geometry
{
    double length[dim];
};

struct Material
{
    double lambda;
    double mu;
    double C[3][3][3][3];
    void init_elasticity_tensor();
};

struct DirichletBC
{
    unsigned int node;
    unsigned int dof;
    double val;
};

template<unsigned int dim>
class Mesh
{
public:
    void generate(Geometry<dim> geom);

    unsigned int Nel;
    unsigned int Nq;
    unsigned int Nen;
    MeshType mesh_type;
    ShapeFcnType shp_func_type;
    QuadratureType quadrature_type;
    unsigned int Nnd;
    unsigned int divs[dim];

    Eigen::MatrixXd Qpts;
    Eigen::VectorXd wts;
    Eigen::MatrixXd N;
    Eigen::MatrixXd dNdxi[dim];

    std::vector<coords> nodal_coords;
    std::vector<connectivity> elem_conn;
};

template<unsigned int dim>
class Mesh2
{
public:
    void generate(Geometry<dim> geom);

    unsigned int Nel;
    unsigned int Nen;
    MeshType mesh_type;
    ShapeFcnType shp_func_type;
    unsigned int Nnd;
    unsigned int divs[dim];

    std::vector<coords> nodal_coords;
    std::vector<connectivity> elem_conn;
};

struct SolverControls
{
    AnalysisType antype;
    unsigned int nsteps;
    unsigned int ninc_max;
    double nr_tol;
};

void print_mat(Eigen::MatrixXd& A, std::string filename);
void print_vec(Eigen::VectorXd& A, std::string filename);
void fem_to_vtk_vector(std::string file, std::vector<connectivity>& ec, std::vector<coords>& nc, Eigen::VectorXd& d);