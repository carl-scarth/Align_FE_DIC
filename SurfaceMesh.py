# Contains useful operations for a FE surface mesh 
import numpy as np

class SurfaceMesh:
    # Contains properties of a finite element surface mesh comprised of quad elements
    def __init__(self, nodes, connectivity, convert_conn = True):
        # nodes = numpy array containing node coordinates
        # connectivity = numpy array containing element connectivity information
        self.nodes = nodes
        self.elements = []
        # Define a QuadElement object for every entry of connectivity
        self.elements = [QuadElement(i, conn_i, convert_conn=convert_conn) for i, conn_i in enumerate(connectivity)]
        
        # get number of nodes and number of elements
        self.n_el = len(self.elements)
        self.n_nodes = self.nodes.shape[0]
        self.is_grid = False # is the mesh a structured grid of elements

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
