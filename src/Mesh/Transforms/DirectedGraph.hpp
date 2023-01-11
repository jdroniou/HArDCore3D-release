// Classes to assist with mesh handling and coarsening

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <assert.h>

#ifndef _DIRECTEDGRAPH_HPP
#define _DIRECTEDGRAPH_HPP

/*!	
 * @defgroup TransformMeshes 
 * @brief Provides tools to transform a mesh (coarsen, change format, etc.)
 */

/*!
 * \addtogroup TransformMeshes
 * @{
 */

// ----------------------------------------------------------------------------
//                            Free functions
// ----------------------------------------------------------------------------
 
/// Type definition for a matrix of unsigned ints
typedef std::vector<std::vector<std::size_t>> Array;

// ----------------------------------------------------------------------------
//                            DirectedFace class definition
// ----------------------------------------------------------------------------

/** The DirectedFace class is the defining element of a DirectedCell. It takes its two end points as parameters and
	contains functions to assist with comparing faces.
 **/

class DirectedFace
{
public:
	/// vector of vertices
    std::vector<std::size_t> V;
	
	/// Null constructor
    DirectedFace();
    
    /// Destructor
    ~DirectedFace();
    
    /**@brief Class Constructor: Initialises the face with the two nodes that form the face
     **/
    DirectedFace(std::vector<std::size_t>);

    /// Add a vertex to face
    void add_vertex(
    	size_t ///< Vertex to be added
    );
};

/// Boolean operation to test if two faces are equal
bool operator==(const DirectedFace &, const DirectedFace &);

// ----------------------------------------------------------------------------
//                            DirectedCell class definition
// ----------------------------------------------------------------------------

/** The DirectedCell class is the defining element of a DirectedGraph. It is defined by a vector of DirectedFace's,
	and a vector of unsigned int's describing the fine cells that make up the cell upon coarsening. The class contains functions
	related to cell manipulation.
 **/

class DirectedCell
{
public:
	/// The faces the cell consists of
    std::vector<DirectedFace> T;
    
    /// The cell ID's of the cell's that have formed this cell
    std::vector<std::size_t> part;
    
    /// Null Constructor
    DirectedCell();
	
	/// Constructor
    DirectedCell(std::size_t);
    
    /// Destructor
    ~DirectedCell();
    
    /**@brief Class Constructor: Initialises the cell with a vector of faces and the cell ID
     **/
    DirectedCell(
    	std::vector<DirectedFace> &, ///< A reference to the faces that the cell consists of
    	std::size_t ///< Initial cell ID
    );
	
	/// Method to append an face to the end of the cell
    void add_face(
    	DirectedFace & ///< Face to be added
    );
	
	/// Remove face at a given position
    void remove_face(
    	std::size_t ///< Position of face to be removed
    );
	
	/** Appends the vector of faces and partition vector of a cell to the
	  *	faces and partition of this cell
	 **/
    void append_cell(
    	DirectedCell & ///< Cell to be added
    );

    bool check_cell();
    
    /// Tests if cell has an face
    bool has_face(
    	DirectedFace ///< Face to check for
    );
	
	/// Prints all the faces of the cell
    void print(); // Useful for debugging
};

// ----------------------------------------------------------------------------
//                            NodeArray class definition
// ----------------------------------------------------------------------------

/** The NodeArray class is generated from the directed graph, and combined with the vertex coordinates
  *	fully describes a mesh.
 **/
 
class CellNodeArray
{
public:
	/// The cell-node array
    Array A;
	
	/// Null constructor
    CellNodeArray();
    
    /// Default constructor
    ~CellNodeArray();
    
    /**@brief Class Constructor: Initialises the NodeArray with the cell-node array
     **/
    CellNodeArray(
    	Array & ///< The cell-node array
    );
    
    /// Test if a given node exists in the cell-node array
    bool node_exists(
    	std::size_t ///< Node to test for
    );
    
    /// Subtracts one from all nodes greater than a given node
    void renum_nodes(
    	std::size_t ///< Node given
    );
    
    /// Prints the cell node array to an out file stream
    void print(
    	std::ofstream *, ///< Pointer to the file stream to print to 
        std::size_t
    );
};

class GraphNodeArray
{
public:
	/// The cell-node array
    std::vector<CellNodeArray> GA;
	
	/// Null constructor
    GraphNodeArray();
    
    /// Default constructor
    ~GraphNodeArray();
    
    /**@brief Class Constructor: Initialises the NodeArray with the cell-node array
     **/
    GraphNodeArray(
    	std::vector<CellNodeArray> & ///< The graph-node array
    );
    
    /// Test if a given node exists in the cell-node array
    bool node_exists(
    	std::size_t ///< Node to test for
    );
    
    /// Subtracts one from all nodes greater than a given node
    void renum_nodes(
    	std::size_t ///< Node given
    );
    
    /// Prints the cell node array to an out file stream
    void print(
    	std::ofstream * ///< Pointer to the file stream to print to 
    );
};

class DirectedGraph
{
public:
	/// Vector of cells that form the graph
    std::vector<DirectedCell> G;
    
    /// Null constructor
    DirectedGraph();
    
    /// Default constructor
    ~DirectedGraph();
    
    /**@brief Class Constructor: Initialises the DirectedGraph with the cells
     **/
    DirectedGraph(
    	std::vector<DirectedCell> & ///< Initial cells
    );
	
	/// Appends a cell to the end of the graph
    void add_cell(
    	DirectedCell & ///< Cell to be added
    );

	/// Removes cell from a given position
    void remove_cell(
    	std::size_t ///< Cell to be removed
    );
	
	/// Test if graph has duplicate faces or unordered cells
    bool test_graph();    
	
	/// Randomise order of cells in graph
    void randomise();
	
	/// Order cells in graph by the number of faces in each cell
    void order();
	
	/// Coarsen the graph by merging cells
    void coarsen();
	
	/// Returns the cell-node array the graph corresponds to
    GraphNodeArray graph_to_array();
    
    /// Returns list of fine cells that each coarse cell consists of
    std::string get_partition();
    
    /// Outputs to a .dat file for gnuplot to read
    void plotfile(
    	std::ofstream *, ///< Pointer to the file stream to print to 
    	std::vector<Eigen::VectorXd> & ///< Vertices
    );

private:
	/* Checks is there exists a shared face, and checks if there exists
	   a shared node at a non-face interface i.e. they touch at a corner.
	   Must be true, false to return true. Also stores a vector of shared 
	   faces to pass to merge_cells - which deletes them.
	 */
	bool can_merge(std::size_t, std::size_t, std::vector<size_t> &);
	
	/* Appends second cell to the first. Removes second cell. Removes 
	   duplicate faces of first cell then orders faces. Handles
	   non-simply connected case.
	 */
    void merge_cells(std::size_t, std::size_t, std::vector<size_t> &);    
    
    // Random boolean generator for coarsen()
    std::knuth_b rand_engine;
    bool random_bool(double);
};

//@}

#endif
