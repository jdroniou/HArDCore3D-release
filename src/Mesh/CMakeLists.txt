# Libraries

FILE(GLOB source_list mesh.* mesh_builder.* mesh_reader.*)
ADD_LIBRARY(mesh ${source_list})

ADD_LIBRARY(DirectedGraph OBJECT Transforms/DirectedGraph.cpp)

# Transform RF mesh in VTU file
ADD_EXECUTABLE(rf-to-vtu Transforms/MeshToVTU.cpp)
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

# Flatten faces in a mesh
ADD_EXECUTABLE(flatten-faces Transforms/FlattenFaces.cpp)
TARGET_LINK_LIBRARIES(flatten-faces mesh ${Boost_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(flatten-faces PUBLIC ${EIGEN3_INCLUDE_DIR})

