# Contains useful operations for a FE surface mesh 
import numpy as np

class SurfaceMesh:
    # Contains properties of a finite element surface mesh comprised of quad elements
    def __init__(self, nodes = [], connectivity = [], convert_conn = True, from_file = False, file_string = [], node_file = [], el_file = []):
        # nodes = numpy array containing node coordinates
        # connectivity = numpy array containing element connectivity information
        
        # Determine if correct information has been provided to initialise mesh object
        if from_file:
            if not (file_string or (node_file and el_file)):
                raise Exception("You specifed \"from_file = True\", please input names for node_file and el_file, or file_string used to indentify both")
        else:
            if not (nodes and connectivity):
                raise Exception("You specified \"from_file = False\", please input nodes and connectivity matrices")
        
        if from_file:
            nodes, connectivity = get_mesh_from_file(file_string, node_file, el_file)

        self.nodes = nodes
        # Define a QuadElement object for every entry of connectivity
        self.elements = [QuadElement(i, conn_i, convert_conn=convert_conn) for i, conn_i in enumerate(connectivity)]
        
        # get number of nodes and number of elements
        self.n_el = len(self.elements)
        self.n_nodes = self.nodes.shape[0]
        self.is_grid = False # is the mesh a structured grid of elements

        # Calculate the element normals and centroids
        self.get_normals()
        self.get_centroids()

    def define_struct_grid(self, n_x, n_y):
        # Specify the mesh arranaged as a structured  grid, with n_x elements in the x direction and
        # n_y elements in the y direction
        self.is_grid = True
        self.n_x = n_x
        self.n_y = n_y

    def get_el_nodes(self):
        # Determine the nodes connected to every element and store in the objects
        [element.get_nodes(self.nodes) for element in self.elements]

    def get_normals(self):
        # Determine unit normal of each element
        [element.get_normal(self.nodes) for element in self.elements]

    def get_centroids(self):
        # Determine centroids of each element
        self.centroids = np.empty((self.n_el,3))
        for i, element in enumerate(self.elements):
            element.get_centroid(self.nodes)
            self.centroids[i,:] = element.centroid

class QuadElement:
    # Contains properties associated with a quadratic element
    def __init__(self, el_id, connectivity, convert_conn = True):
        self.el_id = el_id
        self.connectivity = connectivity
        # Convert node indices to python indexing if necessary
        if convert_conn:
            self.connectivity = self.connectivity - 1

    def get_nodes(self, nodes_mesh):
        # Extract the nodes connected to the element
        # nodes_mesh = numpy array containing coordinates of all nodes in the mesh, indexed by the entries of self.connectivity
        self.nodes = nodes_mesh[self.connectivity]

    def get_normal(self, nodes_mesh):
        # Calculate a surface normal of the element
        # Determine nodes connected to the element, if they aren't already defined
        if not hasattr(self,"nodes"):
            self.get_nodes(nodes_mesh)

        # Determine two unit vectors on the element surface
        a = self.nodes[1,:] - self.nodes[0,:]
        b = self.nodes[3,:] - self.nodes[0,:]
        a = a/np.sqrt(np.sum(a*a))
        b = b/np.sqrt(np.sum(b*b))    
        # Get unit normal by taking the cross product of a and b
        self.n = np.cross(a,b)

    def get_centroid(self, nodes_mesh):
        # Calculate the element centroid
        # Determine nodes connected to the element, if they aren't already defined
        if not hasattr(self,"nodes"):
            self.get_nodes(nodes_mesh)
        self.centroid = np.sum(self.nodes, axis=0)/4.0


def get_mesh_from_file(file_string = [], node_file = [], el_file = []):
    # Load in elements and nodes from a finite element model of the outer surface of the structure
    # over which the DIC data lies
    # file_string = string used to identify both the location of the nodes and elements, appended by file_string_nodes.csv and file_string_elements.csv
    # node_file = name of the file where the nodes are stored, if the labelling differs from the element file
    # el_file = name of the file where elements are stored, if the labelling differs from the element file
    
    # Has the location of the input files been provided correctly?
    if not file_string and not (node_file and el_file):
        raise Exception("Please provide either \"file_string\" with label for both element and node files, or both separate \"node_file\" and \"el_file\"")
    
    if file_string:
        node_file = file_string + "_nodes.csv"
        el_file = file_string + "_elements.csv"

    # Read element connectivities
    connectivity = np.loadtxt(el_file, dtype = int, delimiter = ',', skiprows = 1)
    # read nodal coordinates
    nodes = np.loadtxt(node_file, delimiter = ',', skiprows = 1)
    return(nodes, connectivity)