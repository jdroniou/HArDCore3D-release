# Libraries

FILE(GLOB source_list mesh.* mesh_builder.* mesh_reader.*)
ADD_LIBRARY(mesh ${source_list})
TARGET_LINK_LIBRARIES(mesh Eigen3::Eigen)

ADD_LIBRARY(DirectedGraph OBJECT Transforms/DirectedGraph.cpp)
TARGET_LINK_LIBRARIES(DirectedGraph Eigen3::Eigen)

# Transform RF mesh in VTU file
ADD_EXECUTABLE(rf-to-vtu Transforms/RFToVTU.cpp)
TARGET_LINK_LIBRARIES(rf-to-vtu mesh plot ${Boost_LIBRARIES} Eigen3::Eigen)

# Translate msh file into RF format
ADD_EXECUTABLE(msh-to-rf Transforms/MSHToRF.cpp)
TARGET_LINK_LIBRARIES(msh-to-rf ${Boost_LIBRARIES} Eigen3::Eigen)

# Mesh Coarsen
ADD_EXECUTABLE(mesh-coarsen $<TARGET_OBJECTS:DirectedGraph> Transforms/MeshCoarsen.cpp)
TARGET_LINK_LIBRARIES(mesh-coarsen ${Boost_LIBRARIES} Eigen3::Eigen)

# Check sanity of mesh
ADD_EXECUTABLE(check-mesh Transforms/CheckMesh.cpp)
TARGET_LINK_LIBRARIES(check-mesh mesh plot ${Boost_LIBRARIES} Eigen3::Eigen)

# Make flat faces in a mesh
ADD_EXECUTABLE(make-flat-faces Transforms/MakeFlatFaces.cpp)
TARGET_LINK_LIBRARIES(make-flat-faces mesh ${Boost_LIBRARIES} Eigen3::Eigen)

# Move vertices of a mesh
ADD_EXECUTABLE(move-vertices Transforms/MoveVertices.cpp)
TARGET_LINK_LIBRARIES(move-vertices mesh ${Boost_LIBRARIES} Eigen3::Eigen)

# Write an FVCA10 format mesh from an RF mesh
ADD_EXECUTABLE(rf-to-fvca10format Transforms/RFToFVCA10format.cpp)
TARGET_LINK_LIBRARIES(rf-to-fvca10format mesh ${Boost_LIBRARIES} Eigen3::Eigen)

