#pragma once
#include "FEMUtility.h"

template <unsigned int dim>
class MixedUPElasticityFEM
{
public:
    void generate_mesh();
    void initialize();
    void solve();
    void post_print(std::string file);

    Geometry<dim> geom;
    Material mat;
    Mesh2<dim> mesh;

    QuadratureType quadrature_type;

    unsigned int Nq;
    Eigen::MatrixXd Qpts;
    Eigen::VectorXd wts;
    Eigen::MatrixXd N;
    Eigen::MatrixXd dNdxi[dim];

    ShapeFcnType shp_func_type_p;
    unsigned int Nen_p;
    Eigen::MatrixXd N_p;
    Eigen::MatrixXd dNdxi_p[dim];
    unsigned int offset_p;

    std::vector<DirichletBC> dirBCs;

private:
    static void get_quadrature_pts(unsigned int& Nq, Eigen::MatrixXd& Qpts, Eigen::VectorXd& wts,
                                   MeshType mesh_type, QuadratureType quadrature_type);
    void init_shp_functions();
    void init_shp_functions_p();
    void assemble();
    void compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de);
    void compute_elem_quantities_p(unsigned int iel, Eigen::MatrixXd& xe, connectivity& conn, Eigen::VectorXui& gidx_p, Eigen::MatrixXd& xp, Eigen::VectorXd& pe);
    void compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac,
                                  Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac);
    void compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl);
    void compute_quad_G(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Gl);
    void apply_dirBC();
    void create_dirichlet_map();
    void reform_full_sol();

private:
    Eigen::MatrixXd K;
    Eigen::VectorXd F;
    Eigen::VectorXd d;

    Eigen::MatrixXd Kd;
    Eigen::VectorXd Fd;
    Eigen::VectorXd dd;

    unsigned int gdofs;
    unsigned int ucdofs;

    std::vector<unsigned int> dirmap;
};