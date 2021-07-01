# Libraries
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp)
ADD_LIBRARY(BoundaryConditions OBJECT BoundaryConditions/BoundaryConditions.cpp)
ADD_LIBRARY(hho-general OBJECT HHO-general/HHO3D.cpp)
ADD_LIBRARY(DirectedGraph OBJECT MeshCoarsen/DirectedGraph.cpp)

# Mesh Coarsen
ADD_EXECUTABLE(mesh-coarsen $<TARGET_OBJECTS:DirectedGraph> MeshCoarsen/MeshCoarsen.cpp)
TARGET_LINK_LIBRARIES(mesh-coarsen ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mesh-coarsen PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh quadrature hybridcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh quadrature hybridcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-locvardiff PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Diffusion Advection Reaction 
ADD_EXECUTABLE(hho-diff-advec-reac $<TARGET_OBJECTS:hho-general> $<TARGET_OBJECTS:TestCase>  $<TARGET_OBJECTS:BoundaryConditions> HHO-diff-advec-reac/HHO_DiffAdvecReac.cpp)
TARGET_LINK_LIBRARIES(hho-diff-advec-reac mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diff-advec-reac PUBLIC ${EIGEN3_INCLUDE_DIR})

# DDR magnetostatics
ADD_EXECUTABLE(ddr-magnetostatics DDR-magnetostatics/ddr-magnetostatics.cpp)
TARGET_LINK_LIBRARIES(ddr-magnetostatics mesh quadrature ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-magnetostatics PUBLIC ${EXTERNAL_INCLUDES})

