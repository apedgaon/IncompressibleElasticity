This C++ code is dependent on Eigen 3.4.0 and will need it downloaded.
	 https://eigen.tuxfamily.org/index.php?title=Main_Page

A path to downloaded eigen directory should be provided to the cmake. 

BUILD INTRUCTIONS:

mkdir build
cd build
cmake -DEIGEN_INCLUDE=/path/to/eigen -DCMAKE_BUILD_TYPE=Release ..
cmake --build .



CODE STRUCTURE:
- InCompProblem.cpp - Contains MAIN function and Part 1, 2, 3 problem definitions
- FEMUtility.h, FEMUtility.cpp - Contains Shape fcn, quadrature, elasticty tensor,
				 mesh generation, printing and vtk utility
- ElasticityFEM.h, ElasticityFEM.cpp - Contains Classical Elasticity FEM Class declaration
				       and implementation
- MixedUPElasticityFEM.h, MixedUPElasticityFEM.cpp - Contains Mixed u-P (2-field) Elasticity FEM Class
						     declaration and implementaion
						     (currently can only handle incompressible materials)

