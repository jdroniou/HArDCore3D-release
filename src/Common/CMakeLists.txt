# Library: basis

FILE(GLOB source_list basis.* localdofspace.* globaldofspace.* variabledofspace.*)
ADD_LIBRARY(basis ${source_list})
TARGET_LINK_LIBRARIES(basis quadrature Eigen3::Eigen)

