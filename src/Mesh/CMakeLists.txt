# Libraries

FILE(GLOB source_list mesh.* mesh_builder.* mesh_reader.*)
ADD_LIBRARY(mesh ${source_list})

ADD_LIBRARY(DirectedGraph OBJECT Transforms/DirectedGraph.cpp)

# Transform RF mesh in VTU file
ADD_EXECUTABLE(rf-to-vtu Transforms/RFToVTU.cpp)
TARGET_LINK_LIBRARIES(rf-to-vtu mesh plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(rf-to-vtu PUBLIC ${EIGEN3_INCLUDE_DIR})

# Translate msh file into RF format
ADD_EXECUTABLE(msh-to-rf Transforms/MSHToRF.cpp)
TARGET_LINK_LIBRARIES(msh-to-rf ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(msh-to-rf PUBLIC ${EIGEN3_INCLUDE_DIR})

# Mesh Coarsen
ADD_EXECUTABLE(mesh-coarsen $<TARGET_OBJECTS:DirectedGraph> Transforms/MeshCoarsen.cpp)
TARGET_LINK_LIBRARIES(mesh-coarsen ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(mesh-coarsen PUBLIC ${EIGEN3_INCLUDE_DIR})

# Check sanity of mesh
ADD_EXECUTABLE(check-mesh Transforms/CheckMesh.cpp)
TARGET_LINK_LIBRARIES(check-mesh mesh plot ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(check-mesh PUBLIC ${EIGEN3_INCLUDE_DIR})

# Make flat faces in a mesh
ADD_EXECUTABLE(make-flat-faces Transforms/MakeFlatFaces.cpp)
TARGET_LINK_LIBRARIES(make-flat-faces mesh ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(make-flat-faces PUBLIC ${EIGEN3_INCLUDE_DIR})

# Move vertices of a mesh
ADD_EXECUTABLE(move-vertices Transforms/MoveVertices.cpp)
TARGET_LINK_LIBRARIES(move-vertices mesh ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(move-vertices PUBLIC ${EIGEN3_INCLUDE_DIR})

# Write an FVCA10 format mesh from an RF mesh
ADD_EXECUTABLE(rf-to-fvca10format Transforms/RFToFVCA10format.cpp)
TARGET_LINK_LIBRARIES(rf-to-fvca10format mesh ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(rf-to-fvca10format PUBLIC ${EIGEN3_INCLUDE_DIR})

