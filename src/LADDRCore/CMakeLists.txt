# Library: laddrcore

FILE(GLOB source_list liealgebra.* laddrcore.* laxgrad.* laxcurl.* laxdiv.* lasxgrad.* lasxcurl.* lasxdiv.*)
ADD_LIBRARY(laddrcore ${source_list})
TARGET_LINK_LIBRARIES(laddrcore ddrcore basis mesh quadrature Eigen3::Eigen)
