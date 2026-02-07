#ifndef FRACTUREDOFSPACE_HPP
#define FRACTUREDOFSPACE_HPP

#pragma once

#include<mesh.hpp>
#include<connected_information_generator.hpp>
#include<vsxgrad.hpp>
#include<basis.hpp>
#include <fractureinterpolator.hpp>

namespace HArDCore3D
{
  /*!
   * \addtogroup Common
   * @{
   */

  /// Base class for global DOF spaces. 
  /** The DOFs are organised by increasing geometric entities dimensions: DOFs of vertices, DOFs of edges, DOFs of faces, DOFs of cells.
  This version allows for cases where the number of DOF on each mesh entity varies from one mesh entity to the next */
  
  class FractureDOFSpace : public VariableDOFSpace
  {
  public:
    typedef std::function<Eigen::Vector3d(const Eigen::Vector3d &, const Cell *, const double & lambda)> FunctionType;
    typedef std::function<Eigen::Matrix3d(const Eigen::Vector3d &, const Cell *, const double & lambda)> GradientDisplacementType;
    typedef std::function<Eigen::Matrix3d(const double &, const double &, const Eigen::Vector3d &, const Cell *)> GradientDisplacementType1;
    typedef std::function<double(const Face *, const double & lambda)> FrictionType;

    /// Constructor
    FractureDOFSpace(
             const VSXGrad & vsx,
            //  const SXGrad & sx,
             const ConnectedInformationGenerator & connectedinfo,
             const FractureInterpolator & fracture_interpolator,
             bool use_threads = true
             );


    //------------------------------------------------------------------------------
    // Accessors
    //------------------------------------------------------------------------------
    
    /// Return the mesh
    const Mesh & mesh() const
    {
      return m_vsxgrad.mesh();
    }

    /// Return the polynomial degree
    const size_t & degree() const
    {
      return m_vsxgrad.degree();
    }

    /// Return the underlying sXgrad space
    const SXGrad & sxgrad() const
    {
      return m_vsxgrad.sXgrad();
    }

    /// Return the total degrees of freedom with respect to fracture_faces
    size_t fracture_dofs() const;

    /// Return the total fracture degrees of freedom with respect to a cell
    Eigen::VectorXd cell_fracture_dofs() const;
    
    /// Return the total degrees of freedom on each fracture_faces
    Eigen::VectorXd fracture_face_dof() const
    {
      size_t n = m_connectedinfo.fracture_face_indices_pointer.size();
      Eigen::VectorXd dofs(n);
      for (size_t i = 0; i < n; ++i)
      {
        int face_idx = m_connectedinfo.fracture_face_indices_pointer.at(i)->global_index();
        dofs(i) = m_vsxgrad.numLocalDofsFace(face_idx);
      }
      return dofs;
    }

    // ----------------------------
    // Local face frame definition
    // ----------------------------
    struct Frame {
        Eigen::Vector3d normal;
        Eigen::Vector3d tangent1;
        Eigen::Vector3d tangent2;
    };

    /// Compute orthonormal frame for a given face i.e. normal, tangent1 and tangent2
    Frame get_frame(const Face& F) const {
        Frame frame;
        frame.normal = F.normal();
        frame.tangent1 = F.face_tangent().col(0);
        frame.tangent2 = F.face_tangent().col(1);
        
        return frame;
    }

    struct FrameDerivative {
        Eigen::VectorXd co1;  // coefficients for t1
        Eigen::VectorXd co2;   // coefficients for t2
        Eigen::MatrixXd Df;     // full Jacobian (3 x n)
    };

    FrameDerivative computeDf(
        const Eigen::Vector3d &t1,
        const Eigen::Vector3d &t2,
        double a, double b,
        const Eigen::VectorXd &grad_a,
        const Eigen::VectorXd &grad_b
    ) const{
        int n = grad_a.size();  // dimension of domain
        FrameDerivative result;
        result.co1.resize(n);
        result.co2.resize(n);
        result.Df.resize(3, n);

        // g(x) = a t1 + b t2
        Eigen::Vector3d g = a * t1 + b * t2;
        double r = g.norm();

        // loop over coordinates x_j
        for (int j = 0; j < n; ++j) {
            double da = grad_a(j);   // ∂_j a
            double db = grad_b(j);   // ∂_j b

            // scalar factor (g·t1) and (g·t2)
            double gdot_t1 = g.dot(t1);
            double gdot_t2 = g.dot(t2);

            double factor = a * da + b * db; // g·∂_j g in coefficients

            // coefficients
            double co1 = (da / r) - (gdot_t1 / (r * r * r)) * factor;
            double co2  = (db / r) - (gdot_t2 / (r * r * r)) * factor;

            result.co1(j) = co1;
            result.co2(j) = co2;

            // construct ∂f/∂x_j = alpha * t1 + beta * t2
            result.Df.col(j) = co1 * t1 + co2 * t2;
        }

        return result;
    }

	   

    /// Interpolator of a continuous function
    Eigen::VectorXd interpolate(
          const FunctionType & q, ///< The function to interpolate
          const double & lambda,
          const int doe_cell = -1, ///< The optional degre of cell quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_face = -1, ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          const int doe_edge = -1 ///< The optional degre of edge quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;

    /// Correction term of Interpolator of a continuous function
    Eigen::VectorXd correction_term(
          const FunctionType & q, ///< The function to interpolate
          const double & lambda,
          const Eigen::VectorXd & v,  ///<Existing interpolator vector
          const int doe_face = -1 ///< The optional degre of face quadrature rules to compute the interpolate. If negative, then 2*degree()+3 will be used.
          ) const;
    

    //------------------------------------------------------------------------------
    // Fracture offsets
    //------------------------------------------------------------------------------


    /// Returns the local fracture offset of the face F for vertex V
    inline size_t localfractureOffset(const Face & F, const Vertex & V) const
    {
      size_t nb_dofs = 0;
      for (int iV=0; iV<F.index_vertex(&V); iV++){
        nb_dofs += m_vsxgrad.numLocalDofsVertex(F.vertex(iV)->global_index());
      }
      return nb_dofs;
    }

    /// Returns the local fracture offset of the face F for edge E
    inline size_t localfractureOffset(const Face & F, const Edge & E) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<F.n_vertices(); iV++){
        nb_dofs += m_vsxgrad.numLocalDofsVertex(F.vertex(iV)->global_index());
      }
      for (int iE=0; iE<F.index_edge(&E); iE++){
        nb_dofs += m_vsxgrad.numLocalDofsEdge(F.edge(iE)->global_index());
      }
      return nb_dofs;
    }

    /// Returns the local fracture offset of the face F
    inline size_t localfractureOffset(const Face & F) const
    {
      size_t nb_dofs = 0;
      for (size_t iV=0; iV<F.n_vertices(); iV++){
        nb_dofs += m_vsxgrad.numLocalDofsVertex(F.vertex(iV)->global_index());
      }
      for (size_t iE=0; iE<F.n_edges(); iE++){
        nb_dofs += m_vsxgrad.numLocalDofsEdge(F.edge(iE)->global_index());
      }
      return nb_dofs;
    }

    /// Returns the fracture offset of the face F
    inline size_t globalfractureOffset(const size_t iF) const
    {
      return ( fracture_face_dof().head(iF) ).sum();
    }

    inline size_t globalfractureOffset(const Face & F) const
    {
      return globalfractureOffset(m_connectedinfo.findFaceIndex(F));
    }

    /// Computes the values of the potential reconstruction at the mesh vertices
    std::vector<VectorRd> computeVertexValues(
                  const Eigen::VectorXd & u   ///< DOFs in the discrete space
                  ) const;

    
    //------------------------------------------------------------------------------
    // Restrictions
    //------------------------------------------------------------------------------

    /// Restrict to the edge (including its vertices) of index iE
    Eigen::VectorXd restrictEdge(size_t iE, const Eigen::VectorXd & vh, size_t iT) const;

    /// Restrict to the face (including vertices and edges) of index iF
    Eigen::VectorXd restrictFace(size_t iF, const Eigen::VectorXd & vh, size_t iT) const;

    /// Restrict to the cell (including vertices, edges, faces, elements and lagrange dofs) of index iT
    Eigen::VectorXd restrictCell(size_t iT, const Eigen::VectorXd & vh) const;

    /// Restrict to the cell (including vertices, edges and faces) of index iT
    Eigen::VectorXd restrictCellnew(size_t iT, const Eigen::VectorXd & vh) const;

    /// Restrict to an edge
    inline Eigen::VectorXd restrict(const Edge & E, const Eigen::VectorXd vh, const Cell & T) const
    {
      return restrictEdge(E.global_index(), vh, T.global_index());
    }

    /// Restrict to a face
    inline Eigen::VectorXd restrict(const Face & F, const Eigen::VectorXd vh, const Cell & T) const
    {
      return restrictFace(F.global_index(), vh, T.global_index());
    }

    /// Restrict to a cell
    inline Eigen::VectorXd restrict(const Cell & T, const Eigen::VectorXd vh) const
    {
      return restrictCell(T.global_index(), vh);
    }

    /// Restrict to a cell
    inline Eigen::VectorXd restrictnew(const Cell & T, const Eigen::VectorXd vh) const
    {
      return restrictCellnew(T.global_index(), vh);
    }

    /// evaluating cell potential locally with its corresponding dof from the global vector
    inline Eigen::Vector3d evaluateCellPotential(size_t iT, const Eigen::Vector3d& x_local, 
                                               const Eigen::VectorXd& dofs) const
    {
        Eigen::MatrixXd basis_at_x(3,m_vsxgrad.cellBases(iT).Polykpo3->dimension());
        for (size_t i = 0; i < m_vsxgrad.cellBases(iT).Polykpo3->dimension(); i++)
        {
          basis_at_x.col(i) = m_vsxgrad.cellBases(iT).Polykpo3->function(i,x_local);
        }
        Eigen::MatrixXd point_eval_matrix = basis_at_x * m_vsxgrad.cellPotential(iT);
        
        // Apply to DOFs to get g(x) ∈ ℝ³
        return point_eval_matrix * dofs;
    }

    /// evaluating cell gradient locally with its corresponding dof from the global vector
    inline Eigen::Matrix3d evaluateCellGradient(size_t iT, const Eigen::Vector3d& x_local, 
                                               const Eigen::VectorXd& dofs) const
    {
        std::vector<Eigen::Matrix3d> M(m_vsxgrad.cellBases(iT).Polyk3x3->dimension());
        for (size_t i = 0; i < m_vsxgrad.cellBases(iT).Polyk3x3->dimension(); i++)
        {
          M[i] = m_vsxgrad.cellBases(iT).Polyk3x3->function(i,x_local);
        }
        Eigen::MatrixXd point_eval_matrix = m_vsxgrad.cellGradient(iT) * dofs;
        
        // Apply to DOFs to get g(x) ∈ ℝ³
        Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < m_vsxgrad.cellBases(iT).Polyk3x3->dimension(); ++i) {
            result += point_eval_matrix(i) * M[i];
        }
        return result;
    }

    /// evaluating face potential locally with its corresponding dof from the global vector
    inline Eigen::Vector3d evaluateFacePotential(size_t iF, const Eigen::Vector3d& x_local, 
                                               const Eigen::VectorXd& dofs) const
    {
        Eigen::MatrixXd basis_at_x(3,m_vsxgrad.faceBases(iF).Polykpo3->dimension());
        for (size_t i = 0; i < m_vsxgrad.faceBases(iF).Polykpo3->dimension(); i++)
        {
          basis_at_x.col(i) = m_vsxgrad.faceBases(iF).Polykpo3->function(i,x_local);
        }
        Eigen::MatrixXd point_eval_matrix = basis_at_x * m_vsxgrad.facePotential(iF);
        
        // Apply to DOFs to get g(x) ∈ ℝ³
        return point_eval_matrix * dofs;
    }

    /// Compute the H^1 discrete norm, with full or symmetric gradient
    std::vector<double> computeH1norms(
                          const std::vector<Eigen::VectorXd> & list_dofs, ///< The list of vectors representing the dofs                          
                          const std::string typegrad="full" ///< computes the norm using the "full" or "symmetric" gradient
                          ) const;    

    /// Compute the L^2 norm, of u_exact-Potential(u_h)
    std::vector<double> computeL2normdisp(
      const FunctionType & exact_solution, ///< The exact solution'
      const double & lambda,
      const Eigen::VectorXd & discrete_solution ///< The vectors representing the coefficient of discrete solution                           
      ) const;   

    /// Compute the L^2 norm, of u_exact-Potential(u_h)
    std::vector<double> computeL2normdispgrad(
      const GradientDisplacementType & exact_gradient, ///< The exact solution'
      const double & lambda,
      const Eigen::VectorXd & discrete_solution ///< The vectors representing the coefficient of discrete solution                           
      ) const; 

    /// Compute the L^2 norm, of jump of u_exact-Potential(u_h)
    std::vector<double> computeL2normdispjump(
      const FunctionType & exact_solution, ///< The exact solution'
      const double & lambda,
      const Eigen::VectorXd & discrete_solution ///< The vectors representing the coefficient of discrete solution                           
      ) const; 
    
    /// Compute the L^2 norm of lambda_n-lambda_{n,h}
    std::vector<double> computeL2normtraction(
      const GradientDisplacementType1 & exact_traction,
      const double & lambda,
      const Eigen::VectorXd & discrete_traction
    )const;

    /// Compute the information of the traction and the discrete solution along the fractures
    std::vector<Eigen::VectorXd> computecontact(
      const FrictionType & friction,
      const Eigen::VectorXd & discrete_traction,
      const double & beta,
      const double & lambda
    )const;

    /// Compute the information of the traction and the discrete solution along the fractures
    Eigen::VectorXd computenormaltraction(
      const std::vector<Face *> & fracture_pointer,
      const Eigen::VectorXd & discrete_traction
    )const;

    /// Compute the L^2 norm of the normal and tangential part of reconstruction on each fracture face
    /// It makes sense only for our particular case of unknown solution for specific fractures given by planes
    std::vector<double> reconstruct_norm(   
      const Eigen::VectorXd & discrete_traction
    )const;


    
    //------------------------------------------------------------------------------
    // Extensions
    //------------------------------------------------------------------------------

    /// Extend a face operator to a cell: starting from a matrix opF acting on the DOFs of F (and its edges and vertices),
    /// redistribute the coefficients in a matrix acting on the DOFs viewed from T.
    Eigen::MatrixXd extendOperator(const Cell & T, const Face & F, const Eigen::MatrixXd & opF) const;

    /// Extend an edge operator to a cell
    Eigen::MatrixXd extendOperator(const Cell & T, const Edge & E, const Eigen::MatrixXd & opE) const;

    /// Extend an edge operator to a face
    Eigen::MatrixXd extendOperator(const Face & F, const Edge & E, const Eigen::MatrixXd & opE) const;

    /// Adds a face operator to a cell operator: distributes and adds the contributions from face F to the operator opT on element T (equivalent to opT + extendOperator(T, F, opF)).
    void extendOperator(const Cell & T, const Face & F, Eigen::Ref<Eigen::MatrixXd> opT, const Eigen::MatrixXd & opF) const;

    /// Adds an edge operator to a cell
    void extendOperator(const Cell & T, const Edge & E, Eigen::MatrixXd & opT, const Eigen::MatrixXd & opE) const;

    /// Adds an edge operator to a face
    void extendOperator(const Face & F, const Edge & E, Eigen::MatrixXd & opF, const Eigen::MatrixXd & opE) const;

    /// Takes an inner product prodF on a face F, and adds its contributions to the inner product prodT on the element T (distributes the contributions according to the DOFs as seen from T)
    void addInnerProductContribution(const Cell & T, const Face & F, Eigen::MatrixXd & prodT, const Eigen::MatrixXd & prodF) const;
    
    //------------------------------------------------------------------------------
    // Global DOF indices for an element T
    //------------------------------------------------------------------------------

    /// Returns a vector listing the global DOFs attached to the element T: vertex DOFs, edge DOFs, face DOFs and element DOFs
    std::vector<size_t> globalDOFIndices(const Cell & T) const;

    /// Returns a vector listing the global DOFs attached to the face F: vertex DOFs, edge DOFs, face DOFs
    std::vector<size_t> globalDOFIndices(const Face & F, const Cell & T) const;  

    /// Returns a vector listing the global DOFs attached to the edge E: vertex DOFs, edge DOFs
    std::vector<size_t> globalDOFIndices(const Edge & E, const Cell & T) const; 

    /// Returns a vector listing the global DOFs attached to the vertex V: vertex DOFs
    std::vector<size_t> globalDOFIndices(const Vertex & V, const Cell & T) const;
    
    private:
    const VSXGrad & m_vsxgrad;
    // const SXGrad &m_sxgrad;
    const ConnectedInformationGenerator & m_connectedinfo;
    const FractureInterpolator & m_fracture_interpolator;
    bool m_use_threads;

      
  };

} // namespace HArDCore3D

#endif
