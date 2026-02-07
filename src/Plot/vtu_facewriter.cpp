// Classes and methods to output VTU files for selected faces
//
// Provides method to plot .vtu file from selected faces and face-centered values.

#include "vtu_facewriter.hpp"
#include <memory>
#include <string>
#include <iostream>


namespace HArDCore3D {

// Constructor
FaceVtuWriter::FaceVtuWriter(const std::vector<Face*>& faces)
    : _faces(faces),
			nfaces(faces.size()),
			nvertices(0),
			vtk_type(nfaces),
			f_vertices(nfaces),
			f_offset(nfaces) {
		    
		// Step 1: Collect all unique vertices from the faces
		size_t local_idx = 0;
		for (size_t iF = 0; iF < nfaces; iF++) {
			Face* face = _faces[iF];
			std::vector<Vertex*> face_verts = face->get_vertices();
			
			for (auto* v : face_verts) {
				if (vertex_map.find(v) == vertex_map.end()) {
					vertex_map[v] = local_idx;
					vertices.push_back(v);
					local_idx++;
				}
			}
		}
		nvertices = vertices.size();
		
		// Step 2: Build connectivity and determine VTK types for each face
		vtk_type.assign(nfaces, 0);
		for (size_t iF = 0; iF < nfaces; iF++) {
			Face* face = _faces[iF];
			std::vector<Vertex*> face_verts = face->get_vertices();
			
			// Determine VTK cell type based on number of vertices
			switch(face_verts.size()) {
				case 3: vtk_type[iF] = 5; break;   // Triangle (VTK_TRIANGLE)
				case 4: vtk_type[iF] = 9; break;   // Quad (VTK_QUAD)
				default: vtk_type[iF] = 7; break;  // Polygon (VTK_POLYGON)
			}
			
			// Store local vertex indices for this face
			f_vertices[iF].resize(face_verts.size());
			for (size_t iV = 0; iV < face_verts.size(); iV++) {
				f_vertices[iF][iV] = vertex_map[face_verts[iV]];
			}
			
			// Compute offset (cumulative number of vertices)
			if (iF == 0) {
				f_offset[0] = f_vertices[iF].size();
			} else {
				f_offset[iF] = f_offset[iF-1] + f_vertices[iF].size();
			}
		}
}


bool FaceVtuWriter::write_header(FILE* pFile){
    fprintf (pFile,"%s","<VTKFile type=\"UnstructuredGrid\"");
    fprintf (pFile,"%s"," version=\"0.1\"");
    fprintf (pFile,"%s\n"," byte_order=\"LittleEndian\">");
    fprintf (pFile,"%s\n","   <UnstructuredGrid>");
    fprintf (pFile,"%s","      <Piece  ");
    fprintf (pFile,"%s%zu%s","NumberOfPoints=\"",nvertices,"\"  ");
    fprintf (pFile,"%s%zu%s\n","NumberOfCells=\"",nfaces,"\">");
		
		return true;
}

bool FaceVtuWriter::write_vertices(FILE* pFile){
    // VERTICES
    fprintf (pFile,"%s\n","         <Points>");
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" NumberOfComponents=\"","3","\" format=\"ascii\">");
		for (auto* v : vertices){
	    // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1} z_{i=1} .... x_{i=n} y_{i=n} z_{i=n}) 
			auto vcoords = v->coords();
			fprintf (pFile,"%f %f %f ",vcoords(0),vcoords(1),vcoords(2));
    }
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");
    fprintf (pFile,"%s\n","         </Points>");

		return true;
}

bool FaceVtuWriter::write_solution(FILE* pFile, const std::vector<double> &sol_face, const std::string &name){
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" Name=\"",name.c_str(),"\" NumberOfComponents=\"1\" format=\"ascii\">");
    // write variables (u_{i=0}, u_{i=1}, ... u_{i=n}) i is the face index
    for(size_t iF = 0; iF < nfaces; iF++){ 
			fprintf (pFile,"%f ",sol_face[iF]);
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

		return true;
}

bool FaceVtuWriter::write_solution(FILE* pFile, const Eigen::VectorXd &sol_face, const std::string &name){
    std::vector<double> sol_v(sol_face.data(), sol_face.data() + sol_face.size());
    return write_solution(pFile, sol_v, name); 
}

bool FaceVtuWriter::write_solution(FILE* pFile, const std::vector<VectorRd> &sol_face, const std::string &name){
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" Name=\"",name.c_str(),"\" NumberOfComponents=\"3\" format=\"ascii\">");
    // write variables (ux_{i=0}, uy_{i=0}, uz_{i=0}, ux_{i=1}, ... uz_{i=n}) i is the face index
    for(size_t iF = 0; iF < nfaces; iF++){ 
			fprintf (pFile,"%f %f %f ",sol_face[iF].x(), sol_face[iF].y(), sol_face[iF].z());
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

		return true;
}

bool FaceVtuWriter::write_cells(FILE* pFile){
    // CONNECTIVITY
    fprintf (pFile,"%s\n","         <Cells>");
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
    // Lists all vertices of all faces
		for (size_t iF = 0; iF < nfaces; iF++){
			for (auto v_idx : f_vertices[iF])	{
	 			fprintf (pFile,"%zu ", v_idx);
			}
		}
    fprintf (pFile,"%s\n","");   
    fprintf (pFile,"%s\n","            </DataArray>");

    
    // CELL OFFSETS.
		// Location of where each face finishes in the connectivity list
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
		for (size_t iF = 0; iF < nfaces; iF++){
			fprintf (pFile,"%zu ",f_offset[iF]);
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

    
    // CELL TYPES 
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">");
		for (size_t iF = 0; iF < nfaces; iF++){
			fprintf (pFile,"%d ",vtk_type[iF]);
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

    fprintf (pFile,"%s\n","         </Cells>");

		return true;
}

bool FaceVtuWriter::write_footer(FILE* pFile){
    fprintf (pFile,"%s\n","      </Piece>  ");
    fprintf (pFile,"%s\n","   </UnstructuredGrid>");
    fprintf (pFile,"%s\n","</VTKFile>");

		return true;
}


} // end of namespace HArDCore3D
