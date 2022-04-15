cmake_minimum_required(VERSION 3.10)

project(InCompElasticityFEM VERSION 1.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(HEADER_FILES ElasticityFEM.h MixedUPElasticityFEM.h FEMUtility.h)
add_executable(InCompElasticityFEM InCompProblem.cpp ElasticityFEM.cpp MixedUPElasticityFEM.cpp FEMUtility.cpp)
target_include_directories(InCompElasticityFEM PUBLIC
							"${EIGEN_INCLUDE}")