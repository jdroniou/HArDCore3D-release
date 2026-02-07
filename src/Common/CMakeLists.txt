# Library: basis

FILE(GLOB source_list basis.* localdofspace.* globaldofspace.* variabledofspace.* BubbleDOFSpace.* fracturedofspace.*)
ADD_LIBRARY(basis ${source_list})
TARGET_LINK_LIBRARIES(basis quadrature Eigen3::Eigen)

