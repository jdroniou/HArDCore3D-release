# Library: laddrcore

FILE(GLOB source_list liealgebra.* laxgrad.* laxcurl.* laxdiv.*)
ADD_LIBRARY(laddrcore ${source_list})
TARGET_LINK_LIBRARIES(laddrcore ddrcore basis mesh quadrature)
