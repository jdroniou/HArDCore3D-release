# Library: mesh

FILE(GLOB source_list mesh.* mesh_builder.* cell.* face.* edge.* vertex.*)
ADD_LIBRARY(mesh ${source_list})
TARGET_LINK_LIBRARIES(mesh quadrature)
