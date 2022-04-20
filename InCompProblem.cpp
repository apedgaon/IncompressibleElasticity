#include <iostream>
#include <cmath>
#include "ElasticityFEM.h"
#include "MixedUPElasticityFEM.h"

constexpr double pi = 3.141592653589793;

void Problem1a()
{
    ElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 1;
    fem.mesh.divs[1] = 1;
    fem.generate_mesh();

    fem.quadrature_type = QuadratureType::one_point;

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_1");
}

void Problem1a2()
{
    ElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 1;
    fem.mesh.divs[1] = 1;
    fem.generate_mesh();

    fem.quadrature_type = QuadratureType::one_point;

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_2");
}

void Problem1b()
{
    ElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 4;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    fem.quadrature_type = QuadratureType::one_point;

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
        else if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_3");
}

void Problem1c()
{
    ElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_quadratic;
    fem.mesh.divs[0] = 4;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
        else if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::two_point;
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_4");
}

void Problem2a()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 1;
    fem.mesh.divs[1] = 1;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::one_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::constant;
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_5");
}

void Problem2a2()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 1;
    fem.mesh.divs[1] = 1;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::one_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::constant;
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_6");
}

void Problem2b()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 4;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
        else if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::two_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::constant;
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_7");
}

void Problem2c()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_quadratic;
    fem.mesh.divs[0] = 4;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
        else if (fabs(y - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 1.0) < 1.0e-8 && fabs(y - 1.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.1 * cos(pi / 4);
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.1 * sin(pi / 4);
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::three_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::lagrange_linear;
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_8");
}

void Problem3a()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 3;
    fem.mesh.divs[1] = 2;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 3.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 1;
            temp.val = -0.2;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::one_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::constant;
    fem.initialize();
    fem.solve();
    fem.post_print("rfem_3_1");
}

void Problem3b()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 6;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 3.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 1;
            temp.val = -0.2;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::two_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::constant;
    fem.initialize();
    fem.solve();
    fem.post_print("rfem_3_2");
}


void Problem3c()
{
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    //fem.mat.lambda = 6e6;    // N cm^-2
    //fem.mat.mu = 2e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::quadrilateral;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_quadratic;
    fem.mesh.divs[0] = 6;
    fem.mesh.divs[1] = 4;
    fem.generate_mesh();

    // Store Dirichlet BC
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        double x = fem.mesh.nodal_coords[idx][0];
        double y = fem.mesh.nodal_coords[idx][1];
        if (fabs(x - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }

        if (fabs(x - 3.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 1;
            temp.val = -0.2;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.quadrature_type = QuadratureType::three_point;  // primal u
    fem.shp_func_type_p = ShapeFcnType::lagrange_linear;
    fem.initialize();
    fem.solve();
    fem.post_print("rfem_3_3");
}

int main(int argc, char* argv[])
{
    Problem1a();
    Problem1a2();
    Problem1b();
    Problem1c();

    Problem2a();
    Problem2a2();
    Problem2b();
    Problem2c();

    Problem3a();
    Problem3b();
    Problem3c();
    return 0;
}