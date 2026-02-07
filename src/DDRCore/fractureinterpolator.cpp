#include "fractureinterpolator.hpp"
using namespace HArDCore3D;

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------

FractureInterpolator::FractureInterpolator(
                   const VSXGrad & vsxgrad,
                   const ConnectedInformationGenerator & info_generator,
                   bool use_threads
                   )
  : VariableDOFSpace(
	     vsxgrad.mesh(),
	     [&] {
        Eigen::VectorXi v(vsxgrad.mesh().n_vertices());
        for (size_t i = 0; i < vsxgrad.mesh().n_vertices(); ++i){
            v[i] = info_generator.vertex_final_connected_size.at(vsxgrad.mesh().vertex(i));
        }
        return Eigen::VectorXi((v).matrix());
        }(),
        [&] {
        Eigen::VectorXi e(vsxgrad.mesh().n_edges());
        for (size_t i = 0; i < vsxgrad.mesh().n_edges(); ++i){
            e[i] = info_generator.edge_final_connected_size.at(vsxgrad.mesh().edge(i));
        }
        return Eigen::VectorXi(PolynomialSpaceDimension<Edge>::Poly(vsxgrad.degree()-1) * e.matrix());
        }(),
        [&] {
        Eigen::VectorXi f(vsxgrad.mesh().n_faces());
        for (size_t i = 0; i < vsxgrad.mesh().n_faces(); ++i){
            f[i] = info_generator.face_final_connected_size.at(vsxgrad.mesh().face(i));
        }
        return Eigen::VectorXi(PolynomialSpaceDimension<Face>::Poly(vsxgrad.degree()-1) * f.matrix());
        }(),
	     PolynomialSpaceDimension<Cell>::Poly(vsxgrad.degree()-1) * Eigen::VectorXi::Ones(vsxgrad.mesh().n_cells())),
    m_vsxgrad(vsxgrad),
    m_info_generator(info_generator),
    m_use_threads(use_threads)
    
{
  // Do nothing
}

Eigen::VectorXd FractureInterpolator::interpolate(const ScalarFunctionType & q, const int doe_cell, const int doe_face, const int doe_edge) const
{
  Eigen::VectorXd qh = Eigen::VectorXd::Zero(dimension());

  // Degrees of quadrature rules
  size_t dqr_cell = (doe_cell >= 0 ? doe_cell : 2 * degree() + 3);
  size_t dqr_face = (doe_face >= 0 ? doe_face : 2 * degree() + 3);
  size_t dqr_edge = (doe_edge >= 0 ? doe_edge : 2 * degree() + 3);

  // Interpolate at vertices
  std::function<void(size_t, size_t)> interpolate_vertices
    = [this, &qh, q](size_t start, size_t end)->void
      {
        for (size_t iV = start; iV < end; iV++) {
          const Vertex & V = *(mesh().vertex(iV));
          size_t start_position=globalOffset(V);
          const auto& vert_ptr=mesh().vertex(iV);
          // std::cout<<numLocalDofsVertex(iV)<<"\n";
          for(size_t i=0;i<numLocalDofsVertex(iV);++i){
            const std::vector<Cell *> &current_component = m_info_generator.vertex_final_connected_pointer.at(vert_ptr)[i];
            qh(start_position+i)=q(mesh().vertex(iV)->coords(),current_component[0]);
          }
        } // for iV
      };
  parallel_for(mesh().n_vertices(), interpolate_vertices, m_use_threads);

  if (degree() > 0) {
    
    // Interpolate at edges
    std::function<void(size_t, size_t)> interpolate_edges
      = [this, &qh, q, &dqr_edge](size_t start, size_t end)->void
        {
          for (size_t iE = start; iE < end; iE++) {
            const Edge & E = *(mesh().edge(iE));
            QuadratureRule quad_dqr_E = generate_quadrature_rule(E, dqr_edge);
            auto basis_Pkmo_E_quad = evaluate_quad<Function>::compute(*edgeBases(iE).Polykmo, quad_dqr_E);
            const auto& edge_ptr=mesh().edge(iE);
            for(size_t i=0;i<numLocalDofsEdge(iE);++i){
              const std::vector<Cell *> &current_component = m_info_generator.edge_final_connected_pointer.at(edge_ptr)[i];
              Cell* local_cell_ptr = current_component[0];
              std::function<double(const Eigen::Vector3d &)> q_local = 
                  [q, local_cell_ptr](const Eigen::Vector3d &x) {
                    return q(x, local_cell_ptr);
                  };
              qh.segment(globalOffset(E)+PolynomialSpaceDimension<Edge>::Poly(degree() - 1)*i, PolynomialSpaceDimension<Edge>::Poly(degree() - 1)) 
              = l2_projection(q_local, *edgeBases(iE).Polykmo, quad_dqr_E, basis_Pkmo_E_quad);
            }
          } // for iE
        };
    parallel_for(mesh().n_edges(), interpolate_edges, m_use_threads);

    // Interpolate at faces
    std::function<void(size_t, size_t)> interpolate_faces
      = [this, &qh, q, &dqr_face](size_t start, size_t end)->void
        {
          for (size_t iF = start; iF < end; iF++) {
            const Face & F = *(mesh().face(iF));
            QuadratureRule quad_dqr_F = generate_quadrature_rule(F, dqr_face);
            auto basis_Pkmo_F_quad = evaluate_quad<Function>::compute(*faceBases(iF).Polykmo, quad_dqr_F);
            const auto& face_ptr=mesh().face(iF);
            for(size_t i=0;i<numLocalDofsFace(iF);++i){
              const std::vector<Cell *> &current_component = m_info_generator.face_final_connected_pointer.at(face_ptr)[i];
              Cell* local_cell_ptr = current_component[0];
              std::function<double(const Eigen::Vector3d &)> q_local = 
                  [q, local_cell_ptr](const Eigen::Vector3d &x) {
                    return q(x, local_cell_ptr);
                  };
              qh.segment(globalOffset(F)+PolynomialSpaceDimension<Face>::Poly(degree() - 1)*i, PolynomialSpaceDimension<Face>::Poly(degree() - 1)) 
              = l2_projection(q_local, *faceBases(iF).Polykmo, quad_dqr_F, basis_Pkmo_F_quad);
            }
          } // for iF
        };
    parallel_for(mesh().n_faces(), interpolate_faces, m_use_threads);

    // Interpolate at cells
    std::function<void(size_t, size_t)> interpolate_cells
      = [this, &qh, q, &dqr_cell](size_t start, size_t end)->void
        {
          for (size_t iT = start; iT < end; iT++) {
            const Cell & T = *(mesh().cell(iT));
            const auto& cell_ptr=mesh().cell(iT);
            QuadratureRule quad_dqr_T = generate_quadrature_rule(T, dqr_cell);
            auto basis_Pkmo_T_quad = evaluate_quad<Function>::compute(*cellBases(iT).Polykmo, quad_dqr_T);
            MonomialCellIntegralsType int_mono_2km2 = IntegrateCellMonomials(T, 2*degree()-2);
              std::function<double(const Eigen::Vector3d &)> q_local = 
                  [q, cell_ptr](const Eigen::Vector3d &x) {
                    return q(x,cell_ptr);
                  };
            qh.segment(globalOffset(T), PolynomialSpaceDimension<Cell>::Poly(degree() - 1)) 
              = l2_projection(q_local, *cellBases(iT).Polykmo, quad_dqr_T, basis_Pkmo_T_quad,GramMatrix(T, *cellBases(iT).Polykmo, int_mono_2km2));
          } // for iT
        };
    parallel_for(mesh().n_cells(), interpolate_cells, m_use_threads);
  }
  return qh;
}
