#include "vtu_writer.hpp"
#include <Eigen/Dense>
#include <mesh_builder.hpp>
#include <iostream>
#include <boost/program_options.hpp>

//#include <quadraturerule.hpp>

/*!
 * \addtogroup TransformMeshes
 * @{
 */

constexpr int dimspace = 3;
using namespace HArDCore3D;

//*** Default mesh
const std::string mesh_dir = "../../../meshes/";
std::string default_mesh = mesh_dir + "Cubic-Cells/RF_fmt/gcube_2x2x2";

/// Main executable (CheckMesh) to check the sanity of an RF mesh
int main(const int argc, const char *argv[])
{

  // Program options
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Performs various checks on the mesh")
    ("mesh,m", boost::program_options::value<std::string>(), "Set the mesh");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  // Display the help options
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  // Select the mesh
  std::string mesh_file = (vm.count("mesh") ? vm["mesh"].as<std::string>() : default_mesh);

  std::cout << "Mesh: " << mesh_file << std::endl;

  // Build the mesh    
  HArDCore3D::MeshBuilder meshbuilder = HArDCore3D::MeshBuilder(mesh_file);
  std::unique_ptr<HArDCore3D::Mesh> mesh_ptr = meshbuilder.build_the_mesh();

  // Plot the mesh for exploration
  HArDCore3D::VtuWriter plotdata(mesh_ptr.get());
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(mesh_ptr -> n_vertices());
  std::string filename = mesh_file.substr(mesh_file.find_last_of("/\\") + 1) + ".vtu";
  std::cout << "Writing vtu file to " << filename << std::endl;
  plotdata.write_to_vtu(filename, zero_vec);

  // Check various properties of the mesh
  bool valid = true;
  constexpr double eps=1e-9; // relative tolerance for tests
  
  // File for logs
  std::ofstream logs("logs-checkmesh.txt");
  
  //
  // Number of boundary and internal elements
  std::cout << "Check that nb of internal and boundary entities match total nb of entities" << std::endl;
  logs << "Check that nb of internal and boundary entities match total nb of entities" << std::endl;
  if(mesh_ptr->n_i_vertices() + mesh_ptr->n_b_vertices() != mesh_ptr->n_vertices())
  {
      logs << "**** There are " << mesh_ptr->n_i_vertices() << " internal vertices, " << mesh_ptr->n_b_vertices() << " boundary vertices and " << mesh_ptr->n_vertices() << " total vertices.\n";
      valid = false;
  }

  if(mesh_ptr->n_i_edges() + mesh_ptr->n_b_edges() != mesh_ptr->n_edges())
  {
      logs << "**** There are " << mesh_ptr->n_i_edges() << " internal edges, " << mesh_ptr->n_b_edges() << " boundary edges and " << mesh_ptr->n_edges() << " total edges.\n";
      valid = false;
  }

  if(mesh_ptr->n_i_faces() + mesh_ptr->n_b_faces() != mesh_ptr->n_faces())
  {
      logs << "**** There are " << mesh_ptr->n_i_faces() << " internal faces, " << mesh_ptr->n_b_faces() << " boundary faces and " << mesh_ptr->n_faces() << " total faces.\n";
      valid = false;
  }
  
  if(mesh_ptr->n_i_cells() + mesh_ptr->n_b_cells() != mesh_ptr->n_cells())
  {
      logs << "**** There are " << mesh_ptr->n_i_cells() << " internal cells, " << mesh_ptr->n_b_cells() << " boundary cells and " << mesh_ptr->n_cells() << " total cells.\n";
      valid = false;
  }

  /////
  ///// Orientation of faces
  std::cout << "Check orientation of internal faces" << std::endl;
  logs << "Check orientation of internal faces" << std::endl;
  for(size_t iF = 0; iF < mesh_ptr->n_i_faces(); ++iF)
  {
      const Face * F = mesh_ptr->i_face(iF);
      if (F->n_cells() != 2){
        logs << "**** Internal face " << iF << " of center " << F->center_mass() << " has " << F->n_cells() << " bordering cells." << std::endl;
        valid = false;
      }
      const Cell * T1 = F->cell(0);
      const Cell * T2 = F->cell(1);
      
      if(T1->face_orientation(T1->index_face(F)) * T2->face_orientation(T2->index_face(F)) != -1)
      {
        logs << "**** Face " << F->global_index() << " of center (" << F->center_mass().transpose() << ") has orientations issues" << std::endl;
        logs << "   cell/center/orientation: " << std::endl;
        logs << "       " << T1->global_index() << " / (" << T1->center_mass().transpose() << ") / " << T1->face_orientation(T1->index_face(F)) << std::endl;
        logs << "       " << T2->global_index() << " / (" << T2->center_mass().transpose() << ") / " << T2->face_orientation(T2->index_face(F)) << std::endl << std::endl;
        valid = false;
      }
  }

  ////
  //// Boundary of faces
  std::cout << "Check boundary of faces" << std::endl;
  logs << "Check boundary of faces" << std::endl;
  double max_error_vol_faces = 0.;
  double max_error_int_normal_faces = 0.;
  double min_measure_F = 1e5;
  double max_measure_F = 0.;
  for (size_t iF = 0; iF < mesh_ptr->n_faces(); iF++){
    const Face & F = *mesh_ptr->face(iF);
    min_measure_F = std::min(min_measure_F, F.measure());
    max_measure_F = std::max(max_measure_F, F.measure());
    // Check if:
    //      2 * int_F 1 = int_F div_F(x-x_F) = \sum_E int_E (x-x_F).nFE
    //      sum_E int_E nTE=0
    //      nFE it outside F
    
    double int_xExF_nFE = 0.;
    VectorRd int_nFE = VectorRd::Zero();
    for (size_t iE = 0; iE < F.n_edges(); iE++){
      const Edge & E = *F.edge(iE);
      const VectorRd nFE = F.edge_orientation(iE) * F.edge_normal(iE);
      double xExF_nFE = (E.center_mass()-F.center_mass()).dot(nFE);

      // Volume via integral on the boundary
      int_xExF_nFE += E.measure() * xExF_nFE;

      // Integral of the normal
      int_nFE += E.measure() * nFE;
      
      // Check if normal it outside F
      if ( xExF_nFE < eps ){
        valid = false;
        logs << "**** Face " << iF << " of center (" << F.center_mass().transpose() << ") has outer normal issues at edge " << E.global_index() << " of center " << E.center_mass().transpose() << std::endl;
        logs << "      (xE-xF).nFE = " << xExF_nFE << std::endl;
      }
    }

    bool check_volume = (std::abs(2 * F.measure() - int_xExF_nFE) > F.measure() * eps);
    bool check_int_normal = (int_nFE.norm() > F.diam() * eps);
    double error_vol_F = std::abs((F.measure() - int_xExF_nFE/2)/F.measure());
    double error_int_normal_F = int_nFE.norm();
    max_error_vol_faces = std::max(max_error_vol_faces, error_vol_F);
    max_error_int_normal_faces = std::max(max_error_int_normal_faces, error_int_normal_F);
    if ( check_volume || check_int_normal ){
      valid = false;
      logs << "**** Face " << iF << " of center (" << F.center_mass().transpose() << ") has ";
      if ( !check_int_normal ){
        logs << "AREA issues:" << std::endl;
        logs << "    area: " << F.measure() << ", computed from integral on boundary: " << int_xExF_nFE/2 << " (rel. delta: " << error_vol_F << ")" << std::endl;
      }else if ( !check_volume ){
        logs << "BOUNDARY issues:" << std::endl;
        logs << "    integral of normal: " << int_nFE.transpose() << " (norm: " << error_int_normal_F << ")" << std::endl;                
      }else{
        logs << "AREA and BOUNDARY issues:" << std::endl;
        logs << "    area: " << F.measure() << ", computed from integral on boundary: " << int_xExF_nFE/2 << " (rel. delta: " << error_vol_F << ")" << std::endl;
        logs << "    integral of normal: " << int_nFE.transpose() << " (norm: " << error_int_normal_F << ")" << std::endl;                
      }
      logs << "      Face reg factor: " << std::pow(F.diam(), 2)/F.measure() << std::endl;
//        logs << "      Face vertices: " << std::endl;
//        for (auto V : F.get_vertices()){
//          logs << V->coords().transpose() << std::endl;
//        }
//        const VectorRd t1 = F.edge(0)->tangent();
//        const VectorRd t2 = F.edge_normal(0);
//        const VectorRd n = F.normal();
//        logs << "      Face vertices in coordinates edge(0).tangent() / edge_normal(0) / normal to face: " << std::endl;
//        const Vertex * V0 = F.edge(0)->vertex(0);
//        for (auto V : F.get_vertices()){
//          logs << (V->coords().transpose() - V0->coords().transpose()).dot(t1) << " ";
//          logs << (V->coords().transpose() - V0->coords().transpose()).dot(t2) << " ";
//          logs << (V->coords().transpose() - V0->coords().transpose()).dot(n) << std::endl;
//        }
    }
  }

  ////
  //// Boundary of elements
  std::cout << "Check boundary of elements" << std::endl;
  logs << "Check boundary of elements" << std::endl;
  double max_error_vol_elements = 0.;
  double max_error_int_normal_elements = 0.;
  for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++){
    const Cell & T = *mesh_ptr->cell(iT);
    // Check if:
    //      3 * int_T 1 = int_T div(x-x_T) = \sum_F int_F (x-x_T).nTF
    //      sum_F int_F nTF=0
    //      nTF it outside T
    
    double int_xFxT_nTF = 0.;
    VectorRd int_nTF = VectorRd::Zero();
    for (size_t iF = 0; iF < T.n_faces(); iF++){
      const Face & F = *T.face(iF);
      const VectorRd nTF = T.face_normal(iF);
      double xFxT_nTF = (F.center_mass()-T.center_mass()).dot(nTF);

      // Volume via integral on the boundary
      int_xFxT_nTF += F.measure() * xFxT_nTF;

//        QuadratureRule qr_F_1 = generate_quadrature_rule(F, 1);
//        for (size_t iqn=0; iqn < qr_F_1.size(); iqn++){
//          int_xFxT_nTF += qr_F_1[iqn].w * (qr_F_1[iqn].vector() - T.center_mass()).dot(nTF);
//        }

      // Integral of the normal
      int_nTF += F.measure() * nTF;
      
      // Check if normal it outside T
      if ( xFxT_nTF < eps ){
        valid = false;
        logs << "**** Cell " << iT << " of center (" << T.center_mass().transpose() << ") has outer normal issues at face " << F.global_index() << std::endl;
        logs << "      (xF-xT).nTF = " << xFxT_nTF << std::endl;
        if (F.is_boundary()){
          logs << "      (face is a boundary face)" << std::endl;
        }else{
          logs << "      (cells on each side of face: " << F.cell(0)->global_index() << "/" << F.cell(1)->global_index() << ")" << std::endl;
        }
      }
    }

    bool check_volume = (std::abs(dimspace * T.measure() - int_xFxT_nTF) > T.measure() * eps);
    bool check_int_normal = (int_nTF.norm() > std::pow(T.diam(), 2) * eps);
    double error_vol_T = std::abs((T.measure() - int_xFxT_nTF/dimspace)/T.measure());
    double error_int_normal_T = int_nTF.norm();
    max_error_vol_elements = std::max(max_error_vol_elements, error_vol_T);
    max_error_int_normal_elements = std::max(max_error_int_normal_elements, error_int_normal_T);
    if ( check_volume || check_int_normal ){
      valid = false;
      logs << "**** Cell " << iT << " of center (" << T.center_mass().transpose() << ") has ";
      if ( !check_int_normal ){
        logs << "VOLUME issues:" << std::endl;
        logs << "    volume: " << T.measure() << ", computed from integral on boundary: " << int_xFxT_nTF/dimspace << " (rel. delta: " << error_vol_T << ")" << std::endl;
      }else if ( !check_volume ){
        logs << "BOUNDARY issues:" << std::endl;
        logs << "    integral of normal: " << int_nTF.transpose() << " (norm: " << error_int_normal_T << ")" << std::endl;                
      }else{
        logs << "VOLUME and BOUNDARY issues:" << std::endl;
        logs << "    volume: " << T.measure() << ", computed from integral on boundary: " << int_xFxT_nTF/dimspace << " (rel. delta: " << error_vol_T << ")" << std::endl;
        logs << "    integral of normal: " << int_nTF.transpose() << " (norm: " << error_int_normal_T << ")" << std::endl;                
      }
    }
  }
  logs.close();

  // Conclusion            
  if (valid){
    std::vector<double> reg = mesh_ptr->regularity();
    std::cout << "Mesh is ok. Regularity factors: " << reg[0] << " " << reg[1] << std::endl;
  }else{
    std::cout << "Mesh is NOT ok, see logs-checkmesh.txt for details." << std::endl;
  }
  std::cout << std::endl << "Faces: max rel. error area: " << max_error_vol_faces << ", max error int normal: " << max_error_int_normal_faces << std::endl;
  std::cout << "Elements: max rel. error vol: " << max_error_vol_elements << ", max error int normal: " << max_error_int_normal_elements << std::endl;
    
  return 0;
}
