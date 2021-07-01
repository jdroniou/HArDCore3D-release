// Classes and methods to ouput various file formats
//
// Currently provides:
//  - Method to plot .vtu file from a mesh and node values.
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

#include "vtu_writer.hpp"
#include <memory>
#include <string>
#include <iostream>


namespace HArDCore3D {
// Class
VtuWriter::VtuWriter(const Mesh* mesh)
    : _mesh(mesh),
			ncells(mesh->n_cells()),
			vtk_type(ncells),
			c_vertices(ncells),
			c_offset(ncells),
			c_f_vertices(ncells),
			c_f_offset(ncells) {
		vtk_type.assign(ncells,0);
		for (size_t iC=0; iC<ncells; iC++){
			Cell* cell = _mesh->cell(iC);
		  // create vector that indicates the "vtk" type of each cell.
			// Here, only two: tetra, or generic polyhedral
			switch(cell->n_faces()) {
				case 4: vtk_type[iC]=10; break;			
				default: vtk_type[iC]=42; break;
			}

			// We need the list of all vertices of all cells, in order, with an indication of where the vertices of each cell finish. We create:
			// 	(1) c_vertices[iC] = vector of vertices in cell iC
			// 	(2) coffset[iC] = cumulated number of vertices in the cells from 0 to iC

			c_vertices[iC].resize(cell->n_vertices());
			c_vertices[iC] = cell->get_vertices();
			if (iC==0){
				c_offset[0]=c_vertices[iC].size();
			}else{
				c_offset[iC]=c_offset[iC-1]+c_vertices[iC].size();
			}

			// The "face stream" in vtu file must present itself, for each cell iC:
			//		[nb of faces in iC] [nb vertices in face 1 of iC] [vertex 1 in face 1 of IC] [vertex 2 in face 1 of iC]... [vertex N in face 1 of iC] [nb vertices in face 2 of iC] etc.
			// To do so, we create:
			//	(1) c_f_vertices[iC][ilF] = vector of vertices of face ilF (local numbering) of cell iC
			// 	(2) c_f_offset[iC] = location of end of face stream for cell iC (cumulated number of "1+nb of faces+sum over faces of nb of vertices in the faces" for all cells between 0 and iC)
			if (iC==0){
				c_f_offset[0]=1+cell->n_faces();
			}else{
				c_f_offset[iC]=c_f_offset[iC-1]+1+cell->n_faces();
			}

			c_f_vertices[iC].resize(cell->n_faces());
			for (size_t ilF = 0; ilF < cell->n_faces(); ilF++){
				Face* face = cell->face(ilF);
				c_f_vertices[iC][ilF].resize(face->n_vertices());
				c_f_vertices[iC][ilF] = face->get_vertices();
				c_f_offset[iC] += face->n_vertices();
			}
		}

}

// Method to plot a vtu file
bool VtuWriter::write_to_vtu(const std::string filename, const Eigen::VectorXd &sol_vrtx) {

  FILE* pFile=fopen(filename.c_str(),"w");
	write_header(pFile);
	write_vertices(pFile);
	write_solution(pFile, sol_vrtx);
	write_cells(pFile);
	write_footer(pFile);
	fclose(pFile);
	return true;
}

bool VtuWriter::write_header(FILE* pFile){
    fprintf (pFile,"%s","<VTKFile type=\"UnstructuredGrid\"");
    fprintf (pFile,"%s"," version=\"0.1\"");
    fprintf (pFile,"%s\n"," byte_order=\"LittleEndian\">");
    fprintf (pFile,"%s\n","   <UnstructuredGrid>");
    fprintf (pFile,"%s","      <Piece  ");
    fprintf (pFile,"%s%zu%s","NumberOfPoints=\"",_mesh->n_vertices(),"\"  ");
    fprintf (pFile,"%s%zu%s\n","NumberOfCells=\"",ncells,"\">");
		
		return true;
}

bool VtuWriter::write_vertices(FILE* pFile){
    // VERTICES
    fprintf (pFile,"%s\n","         <Points>");
    fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" NumberOfComponents=\"","3","\" format=\"ascii\">");
		for (auto& v : _mesh->get_vertices()){
	    // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1} z_{i=1} .... x_{i=n} y_{i=n} z_{i=n}) i is the node index 
			auto vcoords = v->coords();
			fprintf (pFile,"%f %f %f ",vcoords(0),vcoords(1),vcoords(2));
    }
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");
    fprintf (pFile,"%s\n","         </Points>");

		return true;
}

bool VtuWriter::write_solution(FILE* pFile, const Eigen::VectorXd &sol_vrtx){
    // SOLUTIONS 
    fprintf (pFile,"%s\n","         <PointData>");
		// in case we have more than one solution, not needed yet
		for(size_t iS=0; iS<1; iS++){
      fprintf (pFile,"%s%s%s\n","            <DataArray type=\"Float32\" Name=\"","solution","\" NumberOfComponents=\"1\" format=\"ascii\">");
      // write variables (u_{i=0}, u_{i=2}, ... u_{i=n}) i is the nodes index
	    for(size_t iV = 0; iV < _mesh->n_vertices(); iV++){ 
				fprintf (pFile,"%f ",sol_vrtx(iV));
			}
      fprintf (pFile,"%s\n","");
      fprintf (pFile,"%s\n","            </DataArray>");
    }
    fprintf (pFile,"%s\n","         </PointData>");

		return true;
}

bool VtuWriter::write_cells(FILE* pFile){
    // CONNECTIVITY
    fprintf (pFile,"%s\n","         <Cells>");
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
    // Lists all vertices of all nodes
		for (size_t iC = 0; iC < ncells; iC++){
			for (auto& v : c_vertices[iC])	{
	 			fprintf (pFile,"%zu ", v->global_index());
			}
		}
    fprintf (pFile,"%s\n","");   
    fprintf (pFile,"%s\n","            </DataArray>");

    
    // CELL OFFSETS.
		// Location of where each cell finishes in the connectivity list
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
		for (size_t iC = 0; iC < ncells; iC++){
			fprintf (pFile,"%zu ",c_offset[iC]);
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

    
    // CELL TYPES 
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">");
		for (size_t iC = 0; iC < ncells; iC++){
			fprintf (pFile,"%d ",vtk_type[iC]);
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

    // FACE STREAM
		// List of all face stream of each cell, where the face stream of cell iC is
		//		[nb of faces in iC] [nb vertices in face 1 of iC] [vertex 1 in face 1 of IC] [vertex 2 in face 1 of iC]... [vertex N in face 1 of iC] [nb vertices in face 2 of iC] etc.
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">");
		for (size_t iC = 0; iC < ncells; iC++){
			// nb of faces in cell iC
			fprintf (pFile,"%zu ",_mesh->cell(iC)->n_faces());
			for (size_t ilF = 0; ilF < c_f_vertices[iC].size(); ilF++){
				// nb of vertices in face ilF of cell iC
				fprintf (pFile,"%zu ",c_f_vertices[iC][ilF].size()); 
				for (size_t iVl = 0; iVl < c_f_vertices[iC][ilF].size(); iVl++){
					// list of vertices of face ilF in cell iC
					fprintf (pFile,"%zu ",c_f_vertices[iC][ilF][iVl]->global_index());
				}
			}
		}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");

    // FACE OFFSETS
		// Positions of the end of each cell iC in the face stream above
    fprintf (pFile,"%s\n","            <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">");
		for (size_t iC = 0; iC < ncells; iC++){
 			fprintf (pFile,"%zu ",c_f_offset[iC]);
	}
    fprintf (pFile,"%s\n","");
    fprintf (pFile,"%s\n","            </DataArray>");


    fprintf (pFile,"%s\n","         </Cells>");

		return true;
}

bool VtuWriter::write_footer(FILE* pFile){
    fprintf (pFile,"%s\n","      </Piece>  ");
    fprintf (pFile,"%s\n","   </UnstructuredGrid>");
    fprintf (pFile,"%s\n","</VTKFile>");

		return true;
}


} // end of namespace HArDCore3D


