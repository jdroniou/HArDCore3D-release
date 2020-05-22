# TestCase
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp TestCase/BoundaryConditions.cpp)
ADD_LIBRARY(hho-general OBJECT HHO-general/HHO3D.cpp)

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diffusion PUBLIC ${EIGEN3_INCLUDE_DIR})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-locvardiff PUBLIC ${EIGEN3_INCLUDE_DIR})

# Magnetostatics solver
ADD_EXECUTABLE(magnetostatics DDRCore/magnetostatics.cpp)  
TARGET_LINK_LIBRARIES(magnetostatics mesh quadrature ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(magnetostatics PUBLIC ${EXTERNAL_INCLUDES})

# HHO Diffusion Advection Reaction 
ADD_EXECUTABLE(hho-diff-advec-reac $<TARGET_OBJECTS:hho-general> $<TARGET_OBJECTS:TestCase>  HHO-diff-advec-reac/HHO_DiffAdvecReac.cpp)
TARGET_LINK_LIBRARIES(hho-diff-advec-reac mesh quadrature hybridcore plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-diff-advec-reac PUBLIC ${EIGEN3_INCLUDE_DIR})

