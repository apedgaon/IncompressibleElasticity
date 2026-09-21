#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "MixedUPElasticityFEM.h"

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::generate_mesh()
{
    mesh.generate(geom);
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::initialize()
{
    // Initialize Shape Functions for displacements and pressures
    get_quadrature_pts(Nq, Qpts, wts, mesh.mesh_type, quadrature_type);
    init_shp_functions();
    init_shp_functions_p();

    // Initilialize global arrays
    gdofs = mesh.Nnd * dim + mesh.Nel * Nen_p;
    offset_p = mesh.Nnd * dim;
    K = Eigen::MatrixXd::Zero(gdofs, gdofs);
    F = Eigen::VectorXd::Zero(gdofs);
    d = Eigen::VectorXd::Zero(gdofs);
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::solve()
{
    // Assemble K, G and F
    assemble();

    // Apply displacement Dirichlet BC
    apply_dirBC();

    // Solve linear system of equations (non positive definite. LU instead of LDLT)
    auto lu = Kd.fullPivLu();
    dd = lu.solve(Fd);

    // Normwise backward error
    Eigen::VectorXd r = Kd * dd - Fd;
    double scale = Kd.norm() * dd.norm() + Fd.norm();
    double error = (scale > 0.0) ? r.norm() / scale : r.norm();

    std::cout << "Finite: " << dd.allFinite() << '\n'
            << "Rank: " << lu.rank() << "/" << Kd.cols() << '\n'
            << "Backward error: " << error << '\n';
    
    // Reform full solution vector
    reform_full_sol();
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::post_print(std::string file)
{
    // Collect only displacement dofs
    Eigen::VectorXd du(mesh.Nnd * dim);
    for (unsigned int idx = 0; idx < mesh.Nnd * dim; ++idx)
        du(idx) = d(idx);

    // Write vtk file with mesh and displacement solutions
    std::cout << "*** Writing VTK file - " << file.c_str() << ".vtu\n";
    fem_to_vtk_vector(file, mesh.elem_conn, mesh.nodal_coords, du);
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::assemble()
{
    // Initialize local arrays and variables for displacments
    unsigned int ldofs = mesh.Nen * dim;
    Eigen::MatrixXd Kl;                      // Local Stiffness
    Eigen::VectorXd Fl;                      // Local Force
    Eigen::MatrixXd xe(dim, mesh.Nen);       // Local Coordinates of nodes
    Eigen::VectorXui gidx(ldofs);            // Global index of local dofs
    Eigen::VectorXd de(mesh.Nen * dim);      // Local nodal displacements
    Eigen::MatrixXd jac(dim, dim);           // Quadrature Isoparamteric Jacobian
    Eigen::MatrixXd jac_inv(dim, dim);       // Quadrature Isoparamteric Jacobian Inverse
    Eigen::MatrixXd dNdX(mesh.Nen, dim);     // Quadrature dN / dX matrices
    double det_jac;                          // Determinant of Isoparametric Jacobian

    // Initialize local arrays and variables for pressure
    unsigned int ldofs_p = Nen_p;
    Eigen::VectorXui gidx_p(ldofs_p);
    Eigen::MatrixXd Gl;
    Eigen::MatrixXd xp(dim, Nen_p);
    Eigen::VectorXd pe(Nen_p);
    unsigned int iel = 0;                     // Element Number

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Initialize local stiffness and force
        Kl = Eigen::MatrixXd::Zero(ldofs, ldofs);
        Fl = Eigen::VectorXd::Zero(ldofs);
        Gl = Eigen::MatrixXd::Zero(ldofs, ldofs_p);

        // Compute Elemental quantities for diplacements and pressure
        compute_elem_quantities(conn, gidx, xe, de);
        compute_elem_quantities_p(iel, xe, conn, gidx_p, xp, pe);


        // Quadrature Loop
        for (unsigned int iq = 0; iq < Nq; ++iq)
        {
            // Compute Quadrature Derivatives for displacements
            compute_quad_derivatives(iq, xe, jac, jac_inv, dNdX, det_jac);

            // Compute and add non-volumetric stiffness to local stiffness matrix Kl
            compute_quad_stiff(iq, dNdX, det_jac, Kl);

            // Compute and add displacement and pressure coupling to local Gl
            compute_quad_G(iq, dNdX, det_jac, Gl);

            // Not impl compute_quad_force()
        }

        // """ K = [Ks  G; G^T  0],   F = [Fs;  0] """
        // Global Assembly Kl -> Ks in K,   Fl -> Fs in F
        for (unsigned int jdx = 0; jdx < ldofs; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
                K(gidx(idx), gidx(jdx)) += Kl(idx, jdx);

            F(gidx(jdx)) += Fl(jdx);
        }

        // Global Assembly Gl -> G in K,    Gl^T -> G^T in K
        for (unsigned int jdx = 0; jdx < ldofs_p; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
            {
                K(gidx(idx), gidx_p(jdx)) += Gl(idx, jdx);
                K(gidx_p(jdx), gidx(idx)) += Gl(idx, jdx);
            }
        }

        ++iel;
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::get_quadrature_pts(unsigned int& Nq, Eigen::MatrixXd& Qpts, Eigen::VectorXd& wts, MeshType mesh_type, QuadratureType quadrature_type)
{
    // Compute and store quadrature points and weights

    if (mesh_type == MeshType::triangular)
    {
        if (quadrature_type == QuadratureType::one_point)
        {
            Nq = 1;
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = 1.0 / 3;
            Qpts(1, 0) = 1.0 / 3;

            wts.resize(Nq);
            wts(0) = 1.0;
        }
        else if (quadrature_type == QuadratureType::two_point)
        {
            Nq = 3;
            Qpts.resize(2, Nq);
            Qpts(0, 0) = 1.0 / 6;
            Qpts(1, 0) = 1.0 / 6;
            Qpts(0, 1) = 2.0 / 3;
            Qpts(1, 1) = 1.0 / 6;
            Qpts(0, 2) = 1.0 / 6;
            Qpts(1, 2) = 2.0 / 3;

            wts.resize(Nq);
            wts(0) = 1.0 / 3;
            wts(1) = 1.0 / 3;
            wts(2) = 1.0 / 3;
        }
    }
    else if (mesh_type == MeshType::quadrilateral)
    {
        if (quadrature_type == QuadratureType::one_point)
        {
            Nq = static_cast<unsigned int>(pow(1, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(1, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(1, 0);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(1, 0) * GsQuad::weights(1, 0);

        }
        else if (quadrature_type == QuadratureType::two_point)
        {
            Nq = static_cast<unsigned int>(pow(2, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 3) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(2, 1);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(1) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0);
            wts(2) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
            wts(3) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
        }
        else if (quadrature_type == QuadratureType::three_point)
        {
            Nq = static_cast<unsigned int>(pow(3, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 2) = GsQuad::quadrature_pt(3, 0);
            Qpts(0, 3) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 4) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 4) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 5) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 5) = GsQuad::quadrature_pt(3, 1);
            Qpts(0, 6) = GsQuad::quadrature_pt(3, 0);
            Qpts(1, 6) = GsQuad::quadrature_pt(3, 2);
            Qpts(0, 7) = GsQuad::quadrature_pt(3, 1);
            Qpts(1, 7) = GsQuad::quadrature_pt(3, 2);
            Qpts(0, 8) = GsQuad::quadrature_pt(3, 2);
            Qpts(1, 8) = GsQuad::quadrature_pt(3, 2);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(3, 0) * GsQuad::weights(3, 0);
            wts(1) = GsQuad::weights(3, 1) * GsQuad::weights(3, 0);
            wts(2) = GsQuad::weights(3, 2) * GsQuad::weights(3, 0);
            wts(3) = GsQuad::weights(3, 0) * GsQuad::weights(3, 1);
            wts(4) = GsQuad::weights(3, 1) * GsQuad::weights(3, 1);
            wts(5) = GsQuad::weights(3, 2) * GsQuad::weights(3, 1);
            wts(6) = GsQuad::weights(3, 0) * GsQuad::weights(3, 2);
            wts(7) = GsQuad::weights(3, 1) * GsQuad::weights(3, 2);
            wts(8) = GsQuad::weights(3, 2) * GsQuad::weights(3, 2);
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::init_shp_functions()
{
    // Compute and store shape functions and their derivatives at quadrature points for displacements
    if (mesh.mesh_type == MeshType::triangular)
    {
        if (mesh.shp_func_type == ShapeFcnType::lagrange_linear)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_tri(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_tri(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdxi_tri(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
    else if (mesh.mesh_type == MeshType::quadrilateral)
    {
        if (mesh.shp_func_type == ShapeFcnType::lagrange_linear)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
        else if (mesh.shp_func_type == ShapeFcnType::lagrange_quadratic)
        {
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = QuadraticLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = QuadraticLagrange::dNdxi_quad(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = QuadraticLagrange::dNdxi_quad(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::init_shp_functions_p()
{
    // Compute and store shape functions and their derivatives at quadrature points for pressure

    if (mesh.mesh_type == MeshType::triangular)
    {
        if (shp_func_type_p == ShapeFcnType::constant)  // one pressure dof at centroid
        {
            Nen_p = 1;
            N_p.resize(Nen_p, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi_p[idx].resize(Nen_p, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nen_p; ++idx)
                {
                    N_p(idx, jdx) = 1.0;
                    dNdxi_p[0](idx, jdx) = 0.0;
                    dNdxi_p[1](idx, jdx) = 0.0;
                }
            }
        }
    }
    else if (mesh.mesh_type == MeshType::quadrilateral)
    {
        if (shp_func_type_p == ShapeFcnType::constant)   // one pressure dof at centroid
        {
            Nen_p = 1;
            N_p.resize(Nen_p, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi_p[idx].resize(Nen_p, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nen_p; ++idx)
                {
                    N_p(idx, jdx) = 1.0;
                    dNdxi_p[0](idx, jdx) = 0.0;
                    dNdxi_p[1](idx, jdx) = 0.0;
                }
            }
        }
        else if (shp_func_type_p == ShapeFcnType::lagrange_linear)  // four pressure dofs at nodes
        {
            Nen_p = 4;
            N_p.resize(Nen_p, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi_p[idx].resize(Nen_p, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nen_p; ++idx)
                {
                    N_p(idx, jdx) = LinearLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi_p[0](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 0, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi_p[1](idx, jdx) = LinearLagrange::dNdxi_quad(idx, 1, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de)
{
    // Computes elemental nodal coordinates, displacements and its global index
    for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            auto loc_idx = dim * ind + idim;
            gidx(loc_idx) = dim * conn[ind] + idim;
            xe(idim, ind) = mesh.nodal_coords[conn[ind]][idim];
            de(loc_idx) = d(gidx(loc_idx));
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::compute_elem_quantities_p(unsigned int iel, Eigen::MatrixXd& xe, connectivity& conn,
                                                          Eigen::VectorXui& gidx_p, Eigen::MatrixXd& xp, Eigen::VectorXd& pe)
{
    // Computes elemental pressure, pressure dof coordinates(not impl) and its global index
    for (unsigned int ip = 0; ip < Nen_p; ++ip)
    {
        gidx_p(ip) = offset_p + iel * Nen_p + ip;
        if (shp_func_type_p == ShapeFcnType::constant)
        {
            //for (unsigned int idim = 0; idim < dim; ++idim)
            //    xp(idim, ip) = 1.0 / 4 * (xe(idim, 0) + xe(idim, 1) + xe(idim, 2) + xe(idim, 3));
        }
        
        pe(ip) = d(gidx_p(ip));
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac, Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac)
{
    // Compute isoparametric Jacobian and Inverse
    jac.setZero();
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                jac(idim, jdim) += xe(idim, ind) * dNdxi[jdim](ind, iq);
        }
    }

    det_jac = jac.determinant();
    jac_inv = jac.inverse();

    // Compute dN/dX
    Eigen::MatrixXd dNdxi_q(mesh.Nen, dim);
    for (unsigned int idim = 0; idim < dim; ++idim)
    {
        for (unsigned int ishp = 0; ishp < mesh.Nen; ++ishp)
            dNdxi_q(ishp, idim) = dNdxi[idim](ishp, iq);
    }

    dNdX = dNdxi_q * jac_inv;
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl)
{
    // Compute Local Deviatoric Stiffness -> integral (sym(w_i,j)  2 mu  sym(u_i,j)) dV
    // """ integral (sym(w_i,j)  2 mu  sym(u_i,j)) dV = integral (w_i,j  mu  u_i,j) dV + integral (w_i,j  mu  u_j,i) dV

    unsigned int loc_idx;
    unsigned int loc_kdx;
    double temp;
    //double mu_rescale = mat.mu;

    // Need this rescaling to avoid nearly zero rows associated with [G^T 0]
    double mu_rescale = 1.0; // K/mu ; p = mu * p_hat

    // integral (w_i,j  mu  u_i,j) dV
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                loc_kdx = dim * bshp + idim;
                loc_idx = dim * ashp + idim;
                temp = 0.0;
                for (unsigned int jdim = 0; jdim < dim; ++jdim)
                    temp += dNdX(ashp, jdim) * mu_rescale *  dNdX(bshp, jdim);

                Kl(loc_idx, loc_kdx) += temp * wts(iq) * det_jac;
            }
        }
    }

    // integral (w_i,j  mu  u_j,i) dV
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int jdim = 0; jdim < dim; ++jdim)
        {
            loc_kdx = dim * bshp + jdim;
            for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
            {
                for (unsigned int idim = 0; idim < dim; ++idim)
                {
                    loc_idx = dim * ashp + idim;
                    Kl(loc_idx, loc_kdx) += dNdX(ashp, jdim) * mu_rescale * dNdX(bshp, idim) * wts(iq) * det_jac;
                }
            }
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::compute_quad_G(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Gl)
{
    // Compute displacement and pressure coupling matrix Gl = integral (w_i,i   p) dV

    unsigned int loc_idx;
    for (unsigned int bshp = 0; bshp < Nen_p; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                loc_idx = dim * ashp + idim;
                Gl(loc_idx, bshp) += dNdX(ashp, idim) * N_p(bshp, iq) * wts(iq) * det_jac;
            }
        }
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::apply_dirBC()
{
    // Applies Dirchlet BC by forming Kd and Fd by removing constrained dofs from K and F

    create_dirichlet_map();

    Kd = Eigen::MatrixXd::Zero(ucdofs, ucdofs);
    Fd = Eigen::VectorXd::Zero(ucdofs);
    dd.resize(ucdofs);
    for (unsigned int jdx = 0; jdx < ucdofs; ++jdx)
    {
        for (unsigned int idx = 0; idx < ucdofs; ++idx)
            Kd(idx, jdx) = K(dirmap[idx], dirmap[jdx]);

        Fd(jdx) = F(dirmap[jdx]);
    }

    for (auto& dirBC : dirBCs)
    {
        unsigned gidx = dirBC.node * dim + dirBC.dof;
        d(gidx) = dirBC.val;
        for (unsigned int idx = 0; idx < ucdofs; ++idx)
            Fd(idx) -= K(dirmap[idx], gidx) * dirBC.val;
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::create_dirichlet_map()
{
    // Creates a local map of unconstrained dofs
    // dirBCs should be sorted for this to work currently
    unsigned int dirBCsize = static_cast<unsigned int>(dirBCs.size());
    ucdofs = gdofs - dirBCsize;       // unconstrained dofs
    dirmap.reserve(ucdofs);
    unsigned int offset = 0;
    auto dir_itr = dirBCs.begin();
    for (unsigned int idx = 0; idx < mesh.Nnd; ++idx)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            if (dir_itr != dirBCs.end())
            {
                if (idx == dir_itr->node && idim == dir_itr->dof)
                {
                    ++dir_itr;
                    continue;
                }
            }

            dirmap.push_back(dim * idx + idim);
        }
    }

    // Needs additional pressure dofs which are by default unconstrained
    for (unsigned int iel = 0; iel < mesh.Nel; ++iel)
    {
        for (unsigned int ien = 0; ien < Nen_p; ++ien)
            dirmap.push_back(offset_p + Nen_p * iel + ien);
    }
}

template<unsigned int dim>
void MixedUPElasticityFEM<dim>::reform_full_sol()
{
    // Constructs full solution vector from solutions at unconstrained dofs
    for (unsigned int idx = 0; idx < ucdofs; ++idx)
        d(dirmap[idx]) = dd(idx);
}

////////////////////////////////////////////////////////////////////////////////
////                            Instantiation                               ////
////////////////////////////////////////////////////////////////////////////////
template class MixedUPElasticityFEM<2>;