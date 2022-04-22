#include <fstream>
#include <iomanip>
#include "FEMUtility.h"

// 1D Element Linear Shape Fcn, Index for correct Shp fcn
double LinearLagrange::N(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = 0.5 * (1 - x);
    else if (index == 1)
        val = 0.5 * (1 + x);

    return val;
}

// 1D Element Derivative of shp fcn
double LinearLagrange::dNdx(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5;
    else if (index == 1)
        val = 0.5;

    return val;
}

// 2D Element (Quad) Linear Shape Fcn, Index for correct Shp fcn
double LinearLagrange::N_quad(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = 0.25 * (1 - xi) * (1 - eta);
    else if (index == 1)
        val = 0.25 * (1 + xi) * (1 - eta);
    else if (index == 2)
        val = 0.25 * (1 + xi) * (1 + eta);
    else if (index == 3)
        val = 0.25 * (1 - xi) * (1 + eta);

    return val;
}


double LinearLagrange::dNdxi_quad(unsigned int index, unsigned int didx, double xi, double eta)
{
    double val = 0.0;
    if (didx == 0) // dNdxi
    {
        if (index == 0)
            val = -0.25 * (1 - eta);
        else if (index == 1)
            val = 0.25 * (1 - eta);
        else if (index == 2)
            val = 0.25 * (1 + eta);
        else if (index == 3)
            val = -0.25 * (1 + eta);
    }
    else if (didx == 1) // dNdeta
    {
        if (index == 0)
            val = -0.25 * (1 - xi);
        else if (index == 1)
            val = -0.25 * (1 + xi);
        else if (index == 2)
            val = 0.25 * (1 + xi);
        else if (index == 3)
            val = 0.25 * (1 - xi);
    }

    return val;
}

// 2D (Tri) Element Linear Shape fcn, index for correct fcn
double LinearLagrange::N_tri(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = 1.0 - xi - eta;
    if (index == 1)
        val = xi;
    else if (index == 2)
        val = eta;

    return val;
}

double LinearLagrange::dNdxi_tri(unsigned int index, unsigned int didx, double xi, double eta)
{
    double val = 0.0;
    if (didx == 0)  // dNdxi
    {
        if (index == 0)
            val = -1.0;
        else if (index == 1)
            val = 1.0;
        else if (index == 2)
            val = 0.0;
    }
    else if (didx == 1)  // dNdeta
    {
        if (index == 0)
            val = -1.0;
        else if (index == 1)
            val = 0.0;
        else if (index == 2)
            val = 1.0;
    }
    return val;
}

// 3D (Hex) Element Linear Shape Fcn, Index for correct Shp fcn
double LinearLagrange::N_hex(unsigned int index, double xi, double eta, double kappa)
{
    double val = 0.0;
    if (index == 0)
        val = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - kappa);
    else if (index == 1)
        val = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - kappa);
    else if (index == 2)
        val = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - kappa);
    else if (index == 3)
        val = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - kappa);
    else if (index == 4)
        val = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + kappa);
    else if (index == 5)
        val = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + kappa);
    else if (index == 6)
        val = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + kappa);
    else if (index == 7)
        val = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + kappa);

    return val;
}

// 3D (Hex) Element Derivatives of Linear Shape Fcn, index for correct Shp fcn, diff_idx for correct derivative
double LinearLagrange::dNdxi_hex(unsigned int index, unsigned int diff_idx, double xi, double eta, double kappa)
{
    double val = 0.0;
    if (diff_idx == 0) // dN/dxi
    {
        if (index == 0)
            val = -0.125 * (1.0 - eta) * (1.0 - kappa);
        else if (index == 1)
            val = 0.125 * (1.0 - eta) * (1.0 - kappa);
        else if (index == 2)
            val = 0.125 * (1.0 + eta) * (1.0 - kappa);
        else if (index == 3)
            val = -0.125 * (1.0 + eta) * (1.0 - kappa);
        else if (index == 4)
            val = -0.125 * (1.0 - eta) * (1.0 + kappa);
        else if (index == 5)
            val = 0.125 * (1.0 - eta) * (1.0 + kappa);
        else if (index == 6)
            val = 0.125 * (1.0 + eta) * (1.0 + kappa);
        else if (index == 7)
            val = -0.125 * (1.0 + eta) * (1.0 + kappa);
    }
    else if (diff_idx == 1) // dN/dEta
    {
        if (index == 0)
            val = -0.125 * (1.0 - xi) * (1.0 - kappa);
        else if (index == 1)
            val = -0.125 * (1.0 + xi) * (1.0 - kappa);
        else if (index == 2)
            val = 0.125 * (1.0 + xi) * (1.0 - kappa);
        else if (index == 3)
            val = 0.125 * (1.0 - xi) * (1.0 - kappa);
        else if (index == 4)
            val = -0.125 * (1.0 - xi) * (1.0 + kappa);
        else if (index == 5)
            val = -0.125 * (1.0 + xi) * (1.0 + kappa);
        else if (index == 6)
            val = 0.125 * (1.0 + xi) * (1.0 + kappa);
        else if (index == 7)
            val = 0.125 * (1.0 - xi) * (1.0 + kappa);
    }
    else if (diff_idx == 2) // dN/dKappa
    {
        if (index == 0)
            val = -0.125 * (1.0 - xi) * (1.0 - eta);
        else if (index == 1)
            val = -0.125 * (1.0 + xi) * (1.0 - eta);
        else if (index == 2)
            val = -0.125 * (1.0 + xi) * (1.0 + eta);
        else if (index == 3)
            val = -0.125 * (1.0 - xi) * (1.0 + eta);
        else if (index == 4)
            val = 0.125 * (1.0 - xi) * (1.0 - eta);
        else if (index == 5)
            val = 0.125 * (1.0 + xi) * (1.0 - eta);
        else if (index == 6)
            val = 0.125 * (1.0 + xi) * (1.0 + eta);
        else if (index == 7)
            val = 0.125 * (1.0 - xi) * (1.0 + eta);
    }

    return val;
}

// 1D Element Quadratic Shp fcn
double QuadraticLagrange::N(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5 * x * (1.0 - x);
    else if (index == 1)
        val = (1.0 - x) * (1.0 + x);
    else if (index == 2)
        val = 0.5 * x * (1.0 + x);

    return val;
}

// 1D Element Derivative of Quadratic Shp fcn
double QuadraticLagrange::dNdx(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5 + x;
    else if (index == 1)
        val = -2.0 * x;
    else if (index == 2)
        val = 0.5 + x;

    return val;
}

double QuadraticLagrange::N_quad(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = 0.25 * xi * (1.0 - xi) * eta * (1.0 - eta);
    else if (index == 1)
        val = -0.25 * xi * (1.0 + xi) * eta * (1.0 - eta);
    else if (index == 2)
        val = 0.25 * xi * (1.0 + xi) * eta * (1.0 + eta);
    else if (index == 3)
        val = -0.25 * xi * (1.0 - xi) * eta * (1.0 + eta);
    else if (index == 4)
        val = -0.5 * (1.0 - xi) * (1.0 + xi) * eta * (1.0 - eta);
    else if (index == 5)
        val = 0.5 * xi * (1.0 + xi) * (1.0 - eta) * (1.0 + eta);
    else if (index == 6)
        val = 0.5 * (1.0 - xi) * (1.0 + xi) * eta * (1.0 + eta);
    else if (index == 7)
        val = -0.5 * xi * (1.0 - xi) * (1.0 - eta) * (1.0 + eta);
    else if (index == 8)
        val = (1.0 - xi) * (1.0 + xi) * (1.0 - eta) * (1.0 + eta);

    return val;
}

double QuadraticLagrange::dNdxi_quad(unsigned int index, unsigned int diff_idx, double xi, double eta)
{
    double val = 0.0;
    if (diff_idx == 0)
    {
        if (index == 0)
            val = 0.25 * (1.0 - 2 * xi) * eta * (1.0 - eta);
        else if (index == 1)
            val = -0.25 * (1.0 + 2 * xi) * eta * (1.0 - eta);
        else if (index == 2)
            val = 0.25 * (1.0 + 2 * xi) * eta * (1.0 + eta);
        else if (index == 3)
            val = -0.25 * (1.0 - 2 * xi) * eta * (1.0 + eta);
        else if (index == 4)
            val = xi * eta * (1.0 - eta);
        else if (index == 5)
            val = 0.5 * (1.0 + 2 * xi) * (1.0 - eta) * (1.0 + eta);
        else if (index == 6)
            val = -xi * eta * (1.0 + eta);
        else if (index == 7)
            val = -0.5 * (1.0 - 2 * xi) * (1.0 - eta) * (1.0 + eta);
        else if (index == 8)
            val = -2.0 * xi * (1.0 - eta) * (1.0 + eta);
    }
    else if (diff_idx == 1)
    {
        if (index == 0)
            val = 0.25 * xi * (1.0 - xi) * (1.0 - 2 * eta);
        else if (index == 1)
            val = -0.25 * xi * (1.0 + xi) * (1.0 - 2 * eta);
        else if (index == 2)
            val = 0.25 * xi * (1.0 + xi) * (1.0 + 2 * eta);
        else if (index == 3)
            val = -0.25 * xi * (1.0 - xi) * (1.0 + 2 * eta);
        else if (index == 4)
            val = -0.5 * (1.0 - xi) * (1.0 + xi) * (1.0 - 2 * eta);
        else if (index == 5)
            val = -xi * (1.0 + xi) * eta;
        else if (index == 6)
            val = 0.5 * (1.0 - xi) * (1.0 + xi) * (1.0 + 2 * eta);
        else if (index == 7)
            val = xi * (1.0 - xi) * eta;
        else if (index == 8)
            val = -2.0 * (1.0 - xi) * (1.0 + xi) * eta;
    }

    return val;
}

// Gauss quadrature points for 1 pt, 2 pt and 3 pt
double GsQuad::quadrature_pt(unsigned int num_pts, unsigned int index)
{
    double retval = 0.0;
    if (num_pts == 1)
        retval = 0.0;
    else if (num_pts == 2)
    {
        if (index == 0)
            retval = -1.0 / sqrt(3.0);
        else if (index == 1)
            retval = 1.0 / sqrt(3.0);
    }
    else if (num_pts == 3)
    {
        if (index == 0)
            retval = -sqrt(3.0 / 5.0);
        else if (index == 1)
            retval = 0.0;
        else if (index == 2)
            retval = sqrt(3.0 / 5.0);
    }

    return retval;
}

// Gauss quadrature weights for 1 pt, 2 pt and 3 pt
double GsQuad::weights(unsigned int num_pts, unsigned int index)
{
    double retval = 0.0;
    if (num_pts == 1)
        retval = 2.0;
    else if (num_pts == 2)
    {
        if (index == 0)
            retval = 1.0;
        else if (index == 1)
            retval = 1.0;
    }
    else if (num_pts == 3)
    {
        if (index == 0)
            retval = 5.0 / 9.0;
        else if (index == 1)
            retval = 8.0 / 9.0;
        else if (index == 2)
            retval = 5.0 / 9.0;
    }

    return retval;
}

// Mesh Generation function (currently implements only hex)
template<unsigned int dim>
inline void Mesh<dim>::generate(Geometry<dim> geom)
{
    switch (quadrature_type)
    {
    // Store quadrature points
    case QuadratureType::two_point:
    {
        if (mesh_type == MeshType::hexahedral)
        {
            Nq = static_cast<unsigned int>(pow(2, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 2) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 3) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 3) = GsQuad::quadrature_pt(2, 0);

            Qpts(0, 4) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 4) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 4) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 5) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 5) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 5) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 7) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 7) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 7) = GsQuad::quadrature_pt(2, 1);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(1) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(2) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1) * GsQuad::weights(2, 0);
            wts(3) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1) * GsQuad::weights(2, 0);

            wts(4) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
            wts(5) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
            wts(6) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
            wts(7) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
        }

    }break;
    }

    switch (shp_func_type)
    {
    case ShapeFcnType::lagrange_linear:
    {
        if (mesh_type == MeshType::hexahedral)
        {
            // Discretization points
            Nen = static_cast<unsigned int>(pow(2, dim));
            unsigned int Nnd_dim[dim];
            Nel = 1;  Nnd = 1;
            double h_dim[dim];
            for (unsigned int idx = 0; idx < dim; ++idx)
            {
                Nel *= divs[idx];
                Nnd_dim[idx] = divs[idx] + 1;
                Nnd *= Nnd_dim[idx];
                h_dim[idx] = geom.length[idx] / divs[idx];
            }

            // Create nodal coordinates
            nodal_coords.reserve(Nnd);
            for (unsigned int kdx = 0; kdx < Nnd_dim[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < Nnd_dim[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < Nnd_dim[0]; ++idx)
                    {
                        coords coordinates(3);
                        coordinates[0] = idx * h_dim[0];
                        coordinates[1] = jdx * h_dim[1];
                        coordinates[2] = kdx * h_dim[2];
                        nodal_coords.push_back(coordinates);
                    }
                }
            }

            // Create element connectivity
            elem_conn.reserve(Nel);
            for (unsigned int kdx = 0; kdx < divs[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < divs[0]; ++idx)
                    {
                        connectivity conn(8);
                        conn[0] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[1] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[2] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[3] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        conn[4] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[5] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[6] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[7] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        elem_conn.push_back(conn);
                    }
                }
            }

            // Store Shape Functions and their derivatives for all quad points
            N.resize(Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_hex(idx, Qpts(0, jdx), Qpts(1, jdx), Qpts(2, jdx));
                    for (unsigned int idim = 0; idim < dim; ++idim)
                        dNdxi[idim](idx, jdx) = LinearLagrange::dNdxi_hex(idx, idim, Qpts(0, jdx), Qpts(1, jdx), Qpts(2, jdx));
                }
            }
        }
        else if (mesh_type == MeshType::triangular)
        {
            Nen = 3;
            Nel = 2 * divs[0] * divs[1];
            unsigned int xnd = divs[0] + 1;
            unsigned int ynd = divs[1] + 1;
            Nnd = xnd * ynd;
            double hx = geom.length[0] / divs[0];
            double hy = geom.length[1] / divs[1];
            nodal_coords.reserve(Nnd);
            for (unsigned int jdx = 0; jdx < ynd; ++jdx)
            {
                for (unsigned int idx = 0; idx < xnd; ++idx)
                {
                    coords coordinates(2);
                    coordinates[0] = idx * hx;
                    coordinates[1] = jdx * hy;
                    nodal_coords.push_back(coordinates);
                }
            }

            elem_conn.reserve(Nel);
            for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
            {
                for (unsigned int idx = 0; idx < divs[0]; ++idx)
                {
                    auto nd0 = jdx * xnd + idx;
                    auto nd1 = jdx * xnd + (idx + 1);
                    auto nd2 = (jdx + 1) * xnd + (idx + 1);
                    auto nd3 = (jdx + 1) * xnd + idx;

                    connectivity conn(3);
                    conn[0] = nd0;
                    conn[1] = nd1;
                    conn[2] = nd3;
                    elem_conn.push_back(conn);

                    conn[0] = nd2;
                    conn[1] = nd3;
                    conn[2] = nd1;
                    elem_conn.push_back(conn);
                }
            }
        }
    }break;

    }
}

template<unsigned int dim>
void Mesh2<dim>::generate(Geometry<dim> geom)
{
    if (mesh_type == MeshType::triangular)
    {
        if (shp_func_type == ShapeFcnType::lagrange_linear)
        {
            Nen = 3;
            Nel = 2 * divs[0] * divs[1];
            unsigned int xnd = divs[0] + 1;
            unsigned int ynd = divs[1] + 1;
            Nnd = xnd * ynd;
            double hx = geom.length[0] / divs[0];
            double hy = geom.length[1] / divs[1];
            nodal_coords.reserve(Nnd);
            for (unsigned int jdx = 0; jdx < ynd; ++jdx)
            {
                for (unsigned int idx = 0; idx < xnd; ++idx)
                {
                    coords coordinates(2);
                    coordinates[0] = idx * hx;
                    coordinates[1] = jdx * hy;
                    nodal_coords.push_back(coordinates);
                }
            }

            elem_conn.reserve(Nel);
            for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
            {
                for (unsigned int idx = 0; idx < divs[0]; ++idx)
                {
                    unsigned int nd0 = jdx * xnd + idx;
                    unsigned int nd1 = jdx * xnd + (idx + 1);
                    unsigned int nd2 = (jdx + 1) * xnd + (idx + 1);
                    unsigned int nd3 = (jdx + 1) * xnd + idx;

                    connectivity conn(3);
                   /* conn[0] = nd0;
                    conn[1] = nd1;
                    conn[2] = nd3;
                    elem_conn.push_back(conn);

                    conn[0] = nd2;
                    conn[1] = nd3;
                    conn[2] = nd1;
                    elem_conn.push_back(conn);*/

                    conn[0] = nd1;
                    conn[1] = nd2;
                    conn[2] = nd0;
                    elem_conn.push_back(conn);

                    conn[0] = nd3;
                    conn[1] = nd0;
                    conn[2] = nd2;
                    elem_conn.push_back(conn);
                }
            }
        }
    }
    else if (mesh_type == MeshType::quadrilateral)
    {
        if (shp_func_type == ShapeFcnType::lagrange_linear)
        {
            // Discretization points
            Nen = static_cast<unsigned int>(pow(2, dim));
            unsigned int Nnd_dim[dim];
            Nel = 1;  Nnd = 1;
            double h_dim[dim];
            for (unsigned int idx = 0; idx < dim; ++idx)
            {
                Nel *= divs[idx];
                Nnd_dim[idx] = divs[idx] + 1;
                Nnd *= Nnd_dim[idx];
                h_dim[idx] = geom.length[idx] / divs[idx];
            }

            // Create nodal coordinates
            nodal_coords.reserve(Nnd);
            for (unsigned int jdx = 0; jdx < Nnd_dim[1]; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nnd_dim[0]; ++idx)
                {
                    coords coordinates(dim);
                    coordinates[0] = idx * h_dim[0];
                    coordinates[1] = jdx * h_dim[1];
                    nodal_coords.push_back(coordinates);
                }
            }

            // Create element connectivity
            elem_conn.reserve(Nel);
            for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
            {
                for (unsigned int idx = 0; idx < divs[0]; ++idx)
                {
                    connectivity conn(4);
                    conn[0] = jdx * Nnd_dim[0] + idx;
                    conn[1] = jdx * Nnd_dim[0] + (idx + 1);
                    conn[2] = (jdx + 1) * Nnd_dim[0] + (idx + 1);
                    conn[3] = (jdx + 1) * Nnd_dim[0] + idx;
                    elem_conn.push_back(conn);
                }
            }

        }
        else if (shp_func_type == ShapeFcnType::lagrange_quadratic)
        {
            // Discretization points
            Nen = static_cast<unsigned int>(pow(3, dim));
            unsigned int Nnd_dim[dim];
            Nel = 1;  Nnd = 1;
            double h_dim[dim];
            for (unsigned int idx = 0; idx < dim; ++idx)
            {
                Nel *= divs[idx];
                Nnd_dim[idx] = 2 * divs[idx] + 1;
                Nnd *= Nnd_dim[idx];
                h_dim[idx] = geom.length[idx] / (2 * divs[idx]);
            }

            // Create nodal coordinates
            nodal_coords.reserve(Nnd);
            coords coordinates(dim);
            for (unsigned int jdx = 0; jdx < Nnd_dim[1]; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nnd_dim[0]; ++idx)
                {
                    coordinates[0] = idx * h_dim[0];
                    coordinates[1] = jdx * h_dim[1];
                    nodal_coords.push_back(coordinates);
                }
            }

            // Create element connectivity
            elem_conn.reserve(Nel);
            connectivity conn(Nen);
            for (unsigned int jdx = 0; jdx < 2 * divs[1]; jdx = jdx + 2)
            {
                for (unsigned int idx = 0; idx < 2 * divs[0]; idx = idx + 2)
                {
                    conn[0] = jdx * Nnd_dim[0] + idx;
                    conn[1] = jdx * Nnd_dim[0] + (idx + 2);
                    conn[2] = (jdx + 2) * Nnd_dim[0] + (idx + 2);
                    conn[3] = (jdx + 2) * Nnd_dim[0] + idx;
                    conn[4] = jdx * Nnd_dim[0] + (idx + 1);
                    conn[5] = (jdx + 1) * Nnd_dim[0] + (idx + 2);
                    conn[6] = (jdx + 2) * Nnd_dim[0] + (idx + 1);
                    conn[7] = (jdx + 1) * Nnd_dim[0] + idx;
                    conn[8] = (jdx + 1) * Nnd_dim[0] + (idx + 1);
                    elem_conn.push_back(conn);
                }
            }
        }
    }
    else if (mesh_type == MeshType::hexahedral)
    {
        if (shp_func_type == ShapeFcnType::lagrange_linear)
        {
            // Discretization points
            Nen = static_cast<unsigned int>(pow(2, dim));
            unsigned int Nnd_dim[dim];
            Nel = 1;  Nnd = 1;
            double h_dim[dim];
            for (unsigned int idx = 0; idx < dim; ++idx)
            {
                Nel *= divs[idx];
                Nnd_dim[idx] = divs[idx] + 1;
                Nnd *= Nnd_dim[idx];
                h_dim[idx] = geom.length[idx] / divs[idx];
            }

            // Create nodal coordinates
            nodal_coords.reserve(Nnd);
            for (unsigned int kdx = 0; kdx < Nnd_dim[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < Nnd_dim[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < Nnd_dim[0]; ++idx)
                    {
                        coords coordinates(3);
                        coordinates[0] = idx * h_dim[0];
                        coordinates[1] = jdx * h_dim[1];
                        coordinates[2] = kdx * h_dim[2];
                        nodal_coords.push_back(coordinates);
                    }
                }
            }

            // Create element connectivity
            elem_conn.reserve(Nel);
            for (unsigned int kdx = 0; kdx < divs[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < divs[0]; ++idx)
                    {
                        connectivity conn(8);
                        conn[0] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[1] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[2] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[3] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        conn[4] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[5] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[6] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[7] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        elem_conn.push_back(conn);
                    }
                }
            }
        }
    }
}

// Forms Elasticty Tensor C_IJKL from isotropic material parameter lambda and mu
void Material::init_elasticity_tensor()
{
    Eigen::Matrix3d del = Eigen::Matrix3d::Identity();
    for (unsigned int idx = 0; idx < 3; ++idx)
    {
        for (unsigned int jdx = 0; jdx < 3; ++jdx)
        {
            for (unsigned int kdx = 0; kdx < 3; ++kdx)
            {
                for (unsigned int ldx = 0; ldx < 3; ++ldx)
                {
                    C[idx][jdx][kdx][ldx] = lambda * del(idx, jdx) * del(kdx, ldx) +
                        mu * (del(idx, kdx) * del(jdx, ldx) + del(idx, ldx) * del(jdx, kdx));
                }
            }
        }
    }
}

// Printing Eigen::Matrix in a specified file
void print_mat(Eigen::MatrixXd& A, std::string filename)
{
    std::ofstream write_stream(filename.c_str());
    if (write_stream.is_open())
    {
        write_stream << std::setprecision(16) << A;
        write_stream.close();
    }
}

void print_vec(Eigen::VectorXd& A, std::string filename)
{
    std::ofstream write_stream(filename.c_str());
    if (write_stream.is_open())
    {
        write_stream << std::setprecision(16) << A;
        write_stream.close();
    }
}

// This function converts solution vector and mesh to vtk file
// This is simply a C++ version of MATLAB code written by Prof. Shiva Rudraraju
void fem_to_vtk_vector(std::string file, std::vector<connectivity>& ec, std::vector<coords>& nc, Eigen::VectorXd& d)
{
    file = file + std::string(".vtu");
    std::ofstream vtk_stream(file.c_str());
    if (vtk_stream.is_open())
    {
        vtk_stream << "<?xml version=\"1.0\"?>\n";
        vtk_stream << "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n";
        vtk_stream << "<UnstructuredGrid>\n";
        auto num_elems = ec.size();
        auto num_nodes = nc.size();
        auto dim_num = nc[0].size();

        // Nodal coords
        vtk_stream << "<Piece  NumberOfPoints=\"" << num_nodes << "\"  NumberOfCells=\"" << num_elems << "\">\n";
        vtk_stream << "<Points>\n";
        vtk_stream << "<DataArray  type=\"Float32\"  NumberOfComponents=\"" << 3 << "\"  format=\"ascii\">\n";
        for (auto& ncoords : nc)
        {
            for (auto ncoord : ncoords)
                vtk_stream << std::setprecision(16) << ncoord << "\t";

            if (dim_num == 2)
                vtk_stream << 0.0;

            vtk_stream << "\n";
        }
        vtk_stream << "</DataArray>\n";

        // Element Connectivity
        vtk_stream << "</Points>\n";
        vtk_stream << "<Cells>\n";
        vtk_stream << "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">\n";
        for (auto& conn : ec)
        {
            for (auto nd_id : conn)
                vtk_stream << nd_id << "\t";

            vtk_stream << "\n";
        }
        vtk_stream << "</DataArray>\n";

        // Offsets
        vtk_stream << "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">\n";
        size_t offset = 0;
        for (auto& conn : ec)
        {
            offset += conn.size();
            vtk_stream << offset << "\n";
        }
        vtk_stream << "</DataArray>\n";

        // Element Types
        vtk_stream << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">\n";
        size_t elem_type = 0;
        for (auto& conn : ec)
        {
            auto elem_order = conn.size();
            if ((dim_num == 2) && (elem_order == 3))
                elem_type = 5;
            else if ((dim_num == 2) && (elem_order == 4))
                elem_type = 9;
            else if ((dim_num == 2) && (elem_order == 9))
                elem_type = 28;
            else if ((dim_num == 3) && (elem_order == 4))
                elem_type = 10;
            else if ((dim_num == 3) && (elem_order == 8))
                elem_type = 12;
            else
                throw std::exception();

            vtk_stream << elem_type << "\n";
        }
        vtk_stream << "</DataArray>\n";
        vtk_stream << "</Cells>\n";

        // vector field
        vtk_stream << "<PointData  Scalars=\"u\">\n";
        vtk_stream << "<DataArray  type=\"Float32\"  Name=\"DisplacementField\"  NumberOfComponents=\"" << 3 << "\"  format=\"ascii\">\n";
        auto field_size = dim_num;
        for (size_t ind = 0; ind < num_nodes; ++ind)
        {
            for (size_t ifld = 0; ifld < field_size; ++ifld)
                vtk_stream << d(ind * field_size + ifld) << "\t";

            if (dim_num == 2)
                vtk_stream << 0.0;

            vtk_stream << "\n";
        }
        vtk_stream << "</DataArray>\n";
        vtk_stream << "</PointData> \n";

        vtk_stream << "</Piece> \n";
        vtk_stream << "</UnstructuredGrid> \n";
        vtk_stream << "</VTKFile> \n";

        vtk_stream.close();
    }
}

/////////////////////////////////////////////////////////////////////////////
//                        Instantiate Mesh<dim>                            //
/////////////////////////////////////////////////////////////////////////////
template class Mesh<2>;
template class Mesh2<2>;
