# Library: ddrcore

FILE(GLOB source_list ddrcore.* xgrad.* xcurl.* xdiv.* serendipity_problem.* sxgrad.* sxcurl.* sxdiv.*)
ADD_LIBRARY(ddrcore ${source_list})
TARGET_LINK_LIBRARIES(ddrcore basis mesh quadrature)
