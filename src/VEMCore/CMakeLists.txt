# Library: vemcore

FILE(GLOB source_list vemcore.* vgrad.* vcurl.* vdiv.*)
ADD_LIBRARY(vemcore ${source_list})
TARGET_LINK_LIBRARIES(vemcore basis mesh quadrature Eigen3::Eigen)
