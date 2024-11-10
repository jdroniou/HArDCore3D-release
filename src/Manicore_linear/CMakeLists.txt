FILE(GLOB cell_list exterior_cell.*)
ADD_LIBRARY(exterior_cell ${cell_list})
TARGET_LINK_LIBRARIES(exterior_cell mesh Eigen3::Eigen)

FILE(GLOB ddr_pec_list ddr_pec.*)
ADD_LIBRARY(ddr_pec ${ddr_pec_list})
TARGET_LINK_LIBRARIES(ddr_pec exterior_cell basis quadrature Eigen3::Eigen)

FILE(GLOB ddr_spaces_list ddr_spaces.*)
ADD_LIBRARY(ddr_spaces ${ddr_spaces_list})
TARGET_LINK_LIBRARIES(ddr_spaces ddr_pec Eigen3::Eigen)

