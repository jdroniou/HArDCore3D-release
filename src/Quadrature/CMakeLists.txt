# Library: quadrature

FILE(GLOB source_list quadraturerule.* triangle_dunavant_rule.* quad3d.* quad3d_face.* quad1d.* tetra_quad.* GMpoly_cell.* GMpoly_face.* GMpoly_edge*)
ADD_LIBRARY(quadrature ${source_list})

