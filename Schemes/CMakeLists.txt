# Libraries
ADD_LIBRARY(TestCase OBJECT TestCase/TestCase.cpp)
TARGET_LINK_LIBRARIES(TestCase Eigen3::Eigen)

ADD_LIBRARY(BoundaryConditions OBJECT BoundaryConditions/BoundaryConditions.cpp)
TARGET_LINK_LIBRARIES(BoundaryConditions Eigen3::Eigen)

ADD_LIBRARY(hho-general OBJECT HHO-general/HHO3D.cpp)
TARGET_LINK_LIBRARIES(hho-general Eigen3::Eigen)

ADD_LIBRARY(MHDTests HHO-MHD/MHDTests.cpp)
TARGET_LINK_LIBRARIES(MHDTests Eigen3::Eigen)

# HHO Diffusion
ADD_EXECUTABLE(hho-diffusion $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-diffusion/HHO_Diffusion.cpp)
TARGET_LINK_LIBRARIES(hho-diffusion mesh hybridcore plot ${EXTERNAL_LIBRARIES})

# HHO Locally variable diffusion
ADD_EXECUTABLE(hho-locvardiff $<TARGET_OBJECTS:TestCase> $<TARGET_OBJECTS:BoundaryConditions> HHO-locvardiff/HHO_LocVarDiff.cpp)
TARGET_LINK_LIBRARIES(hho-locvardiff mesh hybridcore plot ${EXTERNAL_LIBRARIES})

# HHO Diffusion Advection Reaction 
ADD_EXECUTABLE(hho-diff-advec-reac $<TARGET_OBJECTS:hho-general> $<TARGET_OBJECTS:TestCase>  $<TARGET_OBJECTS:BoundaryConditions> HHO-diff-advec-reac/HHO_DiffAdvecReac.cpp)
TARGET_LINK_LIBRARIES(hho-diff-advec-reac mesh hybridcore plot ${Boost_LIBRARIES})

# HHO Full Gradient Diffusion
ADD_EXECUTABLE(hho-fullgradientdiff HHO-full-gradient-diff/hho-fullgradientdiff.cpp)
TARGET_LINK_LIBRARIES(hho-fullgradientdiff hhocore plot BoundaryConditions ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-fullgradientdiff PUBLIC ${EXTERNAL_INCLUDES})

# HHO Brinkman
ADD_EXECUTABLE(hho-brinkman HHO-brinkman/hho-brinkman.cpp)
TARGET_LINK_LIBRARIES(hho-brinkman hhocore plot BoundaryConditions ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(hho-brinkman PUBLIC ${EXTERNAL_INCLUDES})

# DDR magnetostatics (without and with serendipity)
ADD_EXECUTABLE(ddr-magnetostatics DDR-magnetostatics/ddr-magnetostatics.cpp)
TARGET_LINK_LIBRARIES(ddr-magnetostatics ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-magnetostatics PUBLIC ${EXTERNAL_INCLUDES})
ADD_EXECUTABLE(sddr-magnetostatics DDR-magnetostatics/sddr-magnetostatics.cpp)
TARGET_LINK_LIBRARIES(sddr-magnetostatics ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(sddr-magnetostatics PUBLIC ${EXTERNAL_INCLUDES})

# DDR Stokes
ADD_EXECUTABLE(ddr-stokes DDR-stokes/ddr-stokes.cpp)
TARGET_LINK_LIBRARIES(ddr-stokes ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(ddr-stokes PUBLIC ${EXTERNAL_INCLUDES})
ADD_EXECUTABLE(sddr-stokes DDR-stokes/sddr-stokes.cpp)
TARGET_LINK_LIBRARIES(sddr-stokes ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(sddr-stokes PUBLIC ${EXTERNAL_INCLUDES})

# DDR Navier-Stokes
ADD_EXECUTABLE(sddr-navier-stokes $<TARGET_OBJECTS:BoundaryConditions> DDR-navier-stokes/sddr-navier-stokes.cpp)
TARGET_LINK_LIBRARIES(sddr-navier-stokes ddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(sddr-navier-stokes PUBLIC ${EXTERNAL_INCLUDES})

# LADDR Yang-Mills Lagrange multiplier
ADD_EXECUTABLE(laddr-yangmills-lm LADDR-yangmills/laddr-yangmills-lm.cpp)
TARGET_LINK_LIBRARIES(laddr-yangmills-lm laddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(laddr-yangmills-lm PUBLIC ${EXTERNAL_INCLUDES})

# LASDDR Yang-Mills Lagrange multiplier
ADD_EXECUTABLE(lasddr-yangmills LADDR-yangmills/lasddr-yangmills.cpp)
TARGET_LINK_LIBRARIES(lasddr-yangmills laddrcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(lasddr-yangmills PUBLIC ${EXTERNAL_INCLUDES})

# VEM Stokes
ADD_EXECUTABLE(vem-stokes VEM-stokes/vem-stokes.cpp)
TARGET_LINK_LIBRARIES(vem-stokes vemcore plot ${EXTERNAL_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(vem-stokes PUBLIC ${EXTERNAL_INCLUDES})

# HHO MHD
ADD_EXECUTABLE(hho-mhd $<TARGET_OBJECTS:MHDTests> $<TARGET_OBJECTS:BoundaryConditions> HHO-MHD/HHO_MHD.cpp)
TARGET_LINK_LIBRARIES(hho-mhd hybridcore plot mesh ${EXTERNAL_LIBRARIES} Eigen3::Eigen)
TARGET_INCLUDE_DIRECTORIES(hho-mhd PUBLIC ${EXTERNAL_INCLUDES})



