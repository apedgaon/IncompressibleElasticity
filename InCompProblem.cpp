#include <iostream>
#include <cmath>
#include "ElasticityFEM.h"
#include "MixedUPElasticityFEM.h"

// Define pi
constexpr double pi = 3.141592653589793;

// Part 1a: 2 triangle mesh with diagonal displacement (1-field)
void Problem1a()
{
    // Create ElastcityFEM object
    ElasticityFEM<2> fem;

    // Set Geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set Material Parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();   // Constructs C_ijkl

    // Mesh settings: type, shp fcn, divisions
    fem.mesh.mesh_type = MeshType::triangular;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 1;
    fem.mesh.divs[1] = 1;
    fem.generate_mesh(); // generates mesh, nc, ec

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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_1");
}

// Part 1a: 2 triangle mesh with horizontal displacement (1-field)
void Problem1a2()
{
    // Create ElastcityFEM object
    ElasticityFEM<2> fem;

    // Set Geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();   // Consturcts C_ijkl

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_2");
}

// Part 1b: Linear quadrilateral with 1x1 qudrature (1-field)
void Problem1b()
{
    // Create ElastcityFEM object
    ElasticityFEM<2> fem;

    // Set Geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    // Mesh Settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Initialized, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_3");
}

// Part 1c: Quadratic quadrilateral with 2x2 quadrature (1-field)
void Problem1c()
{
    // Create ElasticityFEM object
    ElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor(); // Constructs C_ijkl

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::two_point;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_4");
}

// Part 2a: 2 Triangle mesh with 1x1 quadrature and diagonal displacements (2-field)
void Problem2a()
{
    // Create Mixed UP Elasticity FEM object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::constant;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_5");
}

// Part 2a: 2 Triangle mesh with 1x1 quadrature and horizontal displacements (2-field)
void Problem2a2()
{
    // Create Mixed UP ElasticityFEM object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();  // Constructs C_ijkl

    // Mesh Settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::constant;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_6");
}

// Part 2b: Linear quadrilateral with 2x2 quadrature (2-field)
void Problem2b()
{
    // Create Mixed UP Elasticity FEM object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor(); // Constructs C_ijkl

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::two_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::constant;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_7");
}

// Part 2c: Quadratic Quadrilateral with 3x3 quadrature (2-field)
void Problem2c()
{
    // Create Mixed UP Elastcity FEM Object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 1.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor(); // Constructs C_ijkl

    // Mesh Settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::three_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::lagrange_linear;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("fem_1_8");
}

// Part 3c: Rectangular beam bending with 12 Traingles and 1x1 quadrature (2-field)
void Problem3a()
{
    // Create Mixed UP Elastcity FEM Object
    MixedUPElasticityFEM<2> fem;
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set geometry
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();  // Constructs C_ijkl

    // Mesh Settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::one_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::constant;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("rfem_3_1");
}

// Part 3b: Rectangular beam bending with 24 Linear quadrilaterals and 2x2 quadrature (2-field)
void Problem3b()
{
    // Create Mixed UP Elastcity FEM Object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();   // Construct C_ijkl

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::two_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::constant;

    // Initialize, Solve and Post
    fem.initialize();
    fem.solve();
    fem.post_print("rfem_3_2");
}

// Part 3c: Rectangular beam bending with 12 Quadratic Quadrilaterals and 3x3 quadrature (2-field)
void Problem3c()
{
    // Create Mixed UP Elastcity FEM Object
    MixedUPElasticityFEM<2> fem;

    // Set geometry
    fem.geom.length[0] = 3.0;   // cm
    fem.geom.length[1] = 1.0;   // cm

    // Set material parameters
    fem.mat.lambda = 9.15e8;    // N cm^-2
    fem.mat.mu = 1.83e6;        // N cm^-2
    fem.mat.init_elasticity_tensor();

    // Mesh settings: type, shp fcn, divisions
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

    // Set quadrature type
    fem.quadrature_type = QuadratureType::three_point;

    // Set shape function type for pressure
    fem.shp_func_type_p = ShapeFcnType::lagrange_linear;

    // Initialize, Solve and Post
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