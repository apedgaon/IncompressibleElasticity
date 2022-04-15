#include <iostream>
#include <fstream>
#include <iomanip>
#include "ElasticityFEM.h"

template<unsigned int dim>
void ElasticityFEM<dim>::generate_mesh()
{
    mesh.generate(geom);
}

template<unsigned int dim>
void ElasticityFEM<dim>::initialize()
{
    gdofs = mesh.Nnd * dim;
    K = Eigen::MatrixXd::Zero(gdofs, gdofs);
    F = Eigen::VectorXd::Zero(gdofs);
    d = Eigen::VectorXd::Zero(gdofs);

    init_shp_functions();
}

template<unsigned int dim>
void ElasticityFEM<dim>::solve()
{
    assemble();
    apply_dirBC();
    dd = Kd.ldlt().solve(Fd);
    reform_full_sol();
}

template<unsigned int dim>
void ElasticityFEM<dim>::post_print(std::string file)
{
    //raw_sol_print();

    // Write vtk file with mesh and solutions
    std::cout << "*** Writing VTK file - " << file.c_str() << ".vtu\n";
    fem_to_vtk_vector(file, mesh.elem_conn, mesh.nodal_coords, d);
}


template<unsigned int dim>
void ElasticityFEM<dim>::init_shp_functions()
{
    if (mesh.mesh_type == MeshType::triangular)
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
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_tri(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdeta_tri(idx, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
    }
    if (mesh.mesh_type == MeshType::quadrilateral)
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

        if (mesh.shp_func_type == ShapeFcnType::lagrange_linear)
        {
            // Store Shape Functions and their derivatives for all quad points
            N.resize(mesh.Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(mesh.Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < mesh.Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[0](idx, jdx) = LinearLagrange::dNdxi_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                    dNdxi[1](idx, jdx) = LinearLagrange::dNdeta_quad(idx, Qpts(0, jdx), Qpts(1, jdx));
                }
            }
        }
        else if (mesh.shp_func_type == ShapeFcnType::lagrange_quadratic)
        {
            // Store Shape Functions and their derivatives for all quad points
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
void ElasticityFEM<dim>::assemble()
{
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

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Initialize local stiffness and force
        Kl = Eigen::MatrixXd::Zero(ldofs, ldofs);
        Fl = Eigen::VectorXd::Zero(ldofs);

        // Compute Elemental quantities
        compute_elem_quantities(conn, gidx, xe, de);

        // Quadrature Loop
        for (unsigned int iq = 0; iq < Nq; ++iq)
        {
            compute_quad_derivatives(iq, xe, jac, jac_inv, dNdX, det_jac);

            compute_quad_stiff(iq, dNdX, det_jac, Kl);

            // Not impl compute_quad_force()
        }

        // Global Assembly Kl -> K, Fl -> F
        for (unsigned int jdx = 0; jdx < ldofs; ++jdx)
        {
            for (unsigned int idx = 0; idx < ldofs; ++idx)
                K(gidx(idx), gidx(jdx)) += Kl(idx, jdx);

            F(gidx(jdx)) += Fl(jdx);
        }
    }
}

template<unsigned int dim>
void ElasticityFEM<dim>::compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx,
                                                 Eigen::MatrixXd& xe, Eigen::VectorXd& de)
{
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
void ElasticityFEM<dim>::compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::MatrixXd& jac,
                                                  Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac)
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
void ElasticityFEM<dim>::compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, double det_jac, Eigen::MatrixXd& Kl)
{
    unsigned int loc_idx;
    unsigned int loc_kdx;
    double temp;

    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int kdim = 0; kdim < dim; ++kdim)
        {
            loc_kdx = dim * bshp + kdim;
            for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
            {
                for (unsigned int idim = 0; idim < dim; ++idim)
                {
                    loc_idx = dim * ashp + idim;
                    temp = 0.0;
                    for (unsigned int jdim = 0; jdim < dim; ++jdim)
                    {
                        for (unsigned int ldim = 0; ldim < dim; ++ldim)
                            temp += dNdX(ashp, jdim) * mat.C[idim][jdim][kdim][ldim] * dNdX(bshp, ldim);
                    }

                    Kl(loc_idx, loc_kdx) += temp * wts(iq) * det_jac;
                }
            }
        }
    }
}

template<unsigned int dim>
void ElasticityFEM<dim>::create_dirichlet_map()
{
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
}

template<unsigned int dim>
void ElasticityFEM<dim>::apply_dirBC()
{
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
void ElasticityFEM<dim>::reform_full_sol()
{
    for (unsigned int idx = 0; idx < ucdofs; ++idx)
        d(dirmap[idx]) = dd(idx);
}

template<unsigned int dim>
void ElasticityFEM<dim>::raw_sol_print()
{
    std::ofstream nc("nodal_coords.txt");
    if (nc.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nnd - 1; ++idx)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
                nc << std::setprecision(16) << mesh.nodal_coords[idx][idim] << "\t";

            nc << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
            nc << std::setprecision(16) << mesh.nodal_coords[mesh.Nnd - 1][idim] << "\t";
        nc.close();
    }

    std::ofstream ec("elem_conn.txt");
    if (ec.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nel - 1; ++idx)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                ec << std::setprecision(16) << mesh.elem_conn[idx][ind] << "\t";

            ec << "\n";
        }

        for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
            ec << std::setprecision(16) << mesh.elem_conn[mesh.Nel - 1][ind] << "\t";
        ec.close();
    }

    std::ofstream sol("sol.txt");
    if (sol.is_open())
    {
        for (unsigned int ind = 0; ind < mesh.Nnd - 1; ++ind)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                unsigned int gidx = dim * ind + idim;
                sol << std::setprecision(16) << d(gidx) << "\t";
            }

            sol << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            unsigned int gidx = dim * (mesh.Nnd - 1) + idim;
            sol << std::setprecision(16) << d(gidx) << "\t";
        }

        sol.close();
    }
}

////////////////////////////////////////////////////////////////////////////////
////                            Instantiation                               ////
////////////////////////////////////////////////////////////////////////////////
template class ElasticityFEM<2>;