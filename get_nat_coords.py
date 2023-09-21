# Converts from x, y and z point data to natural coordinates of a finite model
# by projecting points onto the surface and then using a Newton-Raphson solver
# to determine the natural coordinates of the projected points
# Assumes the finite elements are arranged in a regular grid

import numpy as np
from numpy import linalg
from FileSeries import *

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

def get_nat_coords(Files, Mesh, coord_labels = ["x","y","z"]):
    # Read in DIC data from csv files, then apply rotations and translations
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Mesh = Mesh object containing list of nodes and elements, and methods for determining element properties
    # coord_labels = list of strings with labels of the DIC data columns containing the coordinates

    # Calculate the element normals and centroids
    Mesh.get_normals()
    Mesh.get_centroids()

    # Loop through all files in the folder
    for File in Files.files:
        cloud_xyz = File.data[coord_labels].to_numpy() # coordinates        
        # Loop through each of the points in the cloud and project onto the surface of the mesh,
        # then perform Newton-Raphson to determine natural coordinates of this point

        # Consider moving these into another function once the inner workings of the loop have been tidied
        # Maybe something which performs the outer element search?
        xyz_proj = np.empty([File.n_points,3]) # xyz coordinate of projected point
        hr = np.empty([File.n_points,2]) # Array for storing the natural coordinates of each point
        el_ind = [] # list for storing element containing each point
        conv = [] # list containing how many iterations of the Newton-Raphson were required for each point
        for i, point in enumerate(cloud_xyz):
            # Find element with the closest centroid to the point, to use as an intial guess for the element
            min_el = find_closest_centroid(point, Mesh.centroids)
            
            # A for loop with a break might be more sensible... Implement on another iterations

            # Search to determine if the point is found to lie in the chosen element. If  correct element has been 
            # chosen the natural coordinates will lie between +/-1.0, if not the converged values of these coordinatess
            # indicate a direction to move in the grid until the correct element is found, or it is determined that the 
            # point does not lie over the mesh surface
            el_found = False # Boolean to check whether the element the point belongs to has been found
            hr_n_old = np.array([[np.nan], [np.nan]]) # for storing previous converged value of hr when searching
            # If implementing as a for with conditional, rather than while, do this on the first iteration 
            # for an out-of-bounds element to save doing this for every point...
            while not el_found:                
                print(Mesh.elements[min_el].nodes)
                # Project point onto the surface of the current element, and store result
                point_proj = proj_point_on_element(point, Mesh.elements[min_el].centroid, Mesh.elements[min_el].n)
                print(point_proj)
                xyz_proj[i,:] = point_proj
                # Determine natural coordinates for the  point using Newton-Raphson
                # Continue to update function
                (gh, i_conv) = newton_raphson(point_proj, Mesh.elements[min_el].nodes)
                print(gh)
                print(i_conv)
                # replace hr_n with gh below, reflect the change in nomenclature
                # replace j with i_conv below
                # BE CAREFUL!! REMEMBER THAT I ADDED THE DROPNA OPTION TO THE EXAMPLE IN EXPORT 2 TO ENABLE KEEPING NAS TO 
                # HELP WITH INTERPOLATING POINTS ACROSS LOAD STEPS - THIS INVOLVES SOME ADDITION TO THE OUTER LOOP ABOVE
                # AND WILL REQUIRE PASSING A KEYWORD ARGUMENT
                # DO ON NEXT ITERATION
                sadsad

                # Check if the converged natural coordinates are within the element bounds (these
                # should be +/-1.) If not, then another element is selected, and another iteration
                # of the outer while loop is performed
                if np.all((hr_n >= -1) & (hr_n <= 1)):
                    # Condition for exiting the while loop
                    el_found = True
                else:
                    # Job to do on next tidy: replace the while loop with a for loop with a break.
                    # Maybe set loop limit to 5??
                    # There isn't a built in function which will redo the for loop iteration if necessary. 
                    #  It might be tidier to do this using a while loop at the outer level, rather than a for, and
                    #  controlling the incrementation. A break statement inside a conditional will prevent moving to the next iteration. Consider playing with this 
                    # and see which is neater. See:
                    # https://stackoverflow.com/questions/492860/python-restarting-a-loop
                    # https://stackoverflow.com/questions/36573486/redo-for-loop-iteration-in-python

                    # Determine position of current element in the grid to find the next element
                    row_ind = min_el//n_x # Row index. Divide by number of columns and round down to nearest integer
                    col_ind = min_el%n_x # Remainder after the above integer division gives the column index
         
                    # Move along to the next element in the grid based upon the converged value of the 
                    # natural coordinate gives a clue of which direction in the grid the actual element lies in
                    # I've corrected this to match 1,2,3,4 now. Be careful though as this will depend on in which
                    # direction the h and r (g and h) increase with increasing element numbering. I think what I
                    # have is a reasonable default, but it might be worth having this as an input? Perhaps orientation
                    # of the grid relative to g and h is also something to think about
                    if hr_n[0,0] < -1:
                        row_ind = row_ind - 1
                    elif hr_n[0,0] > 1:
                        row_ind = row_ind + 1
                    if hr_n[1,0] < -1:
                        col_ind = col_ind - 1
                    elif hr_n[1,0] > 1:
                        col_ind = col_ind + 1
            
                    # Check to see if either the row or column index are outside the bounds of the grid. 
                    # If so, the point doesn't belong to any element and should be deleted
                    if ((row_ind < 0) | (row_ind >= n_y) | (col_ind < 0) | (col_ind >= n_x)):
                        # Return nans for everything so the points may be deleted outside of the loop
                        # min_el = np.nan - deleted to provent later issues with pandas - just get rid of probably
                        # j = np.nan
                        hr_n = np.array([[np.nan], [np.nan]])
                        el_found = True # Break out of the while loop
                    elif np.any(((hr_n < -1) & (hr_n_old > 1)) | ((hr_n > 1) & (hr_n_old < -1))):
                    # Otherwise, check that the out-of-bounds natural coordinate hasn't flipped signs
                    # from a previous iteration. A flip from hr_n > 1 to  hr_n < -1 or vice-versa  
                    # indicates the point is between the projection of the element face onto the plane
                    # of the point, which can occur if abive a convex curve. If this happens round to 
                    # the boundary of the current element (it doesn't matter which is chosen)
                        if ((hr_n[0,0] < -1) | (hr_n[0,0] > 1)):
                            hr_n[0,0] = round(hr_n[0,0])
                        if ((hr_n[1,0] < -1) | (hr_n[1,0] > 1)):
                            hr_n[1,0] = round(hr_n[1,0])
                        # Is there a neater implementation of the above conditional?? ^^^^^
                        el_found = True # Break out of the while loop
                    else:
                        # If neither, move to the next element in the grid in the direction of the 
                        # converged natural coordinate value
                        min_el = row_ind*n_x + col_ind
                        hr_n_old = hr_n # Store the previous natural coordinate values to track the element search
        
            # Store the index of the element in which the current point sits
            el_ind.append(min_el)
            conv.append(j)
            hr[i,:] = hr_n.squeeze()

        # Remove points which were found to lie outside of the mesh, identified using nans in h and r

    # Dictionary of new column names and their positions in the array and output dataframe respectively
    new_cols = {"x_proj": [0,4], "y_proj": [1,9],"z_proj":[2,14]}
    # Add new entries to the cloud_data dataframe
    hr = pd.DataFrame(hr, columns=["h","r"])
    # Be careful here, if remove na hasn't happened and the row indices don't line up this will add nans!!!
    [cloud_data.insert(loc=value[1], column = key, value = pd.Series(xyz_proj[:,value[0]])) for key, value in new_cols.items()]
    cloud_data = pd.concat([cloud_data, hr],axis=1)
    cloud_data["Element"] = pd.Series(el_ind).astype(int)
    cloud_data["Conv_Iteration"] = pd.Series(conv).astype(int)
    # Identify rows with na value for h and retain for output, dropping h, r, Element and Conv_Teration, as these 
    # do not apply to deleted data
    cloud_data_del = cloud_data[pd.isna(cloud_data["h"])].drop(labels = ["h","r","Element", "Conv_Iteration"],axis=1)
    # Drop rows with nas from main dataframe
    cloud_data = cloud_data.dropna()

    # Switched to 1,2,3,4 rather than 4,1,5,8 so will need to do sanity checks!!!!
    
    # Old code, retain for now
    #nan_ind = np.isnan(hr[:,0])
    #hr = hr[np.logical_not(nan_ind),:]
    #xyz_del = xyz_proj[nan_ind,:]
    #xyz_proj = xyz_proj[np.logical_not(nan_ind),:]
    #cloud_data_del = cloud_data[nan_ind]
    #cloud_data = cloud_data[~nan_ind]
    #el_ind = [el for el in el_ind if not np.isnan(el)]
    #conv = [n for n in conv if not np.isnan(n)] 
    # write adjusted point cloud, and deleted points to csv files
    # [cloud_data.insert(loc=value[1], column = key, value = pd.Series(xyz_proj[:,value[0]])) for key, value in new_cols.items()]
    # [cloud_data_del.insert(loc=value[1], column = key, value = pd.Series(xyz_del[:,value[0]])) for key, value in new_cols.items()]
    # Also append converged, natural coordinates, element index and  number of iterations
    # cloud_data = pd.concat([cloud_data, pd.DataFrame(hr, columns=["h", "r"])],axis=1,ignore_index=True)
    #cloud_data["Element"] = el_ind
    #cloud_data["Convergence Iteration"] = conv
    # WHEN WRITING HERE - HAVE OPTION TO OUTPUT ELEMENT INDICES IN ABAQUS, RATHER THAN PYTHON NOTATION
    cloud_data.to_csv(os.path.join(out_path, filename), sep=",", index=False)
    cloud_data_del.to_csv(os.path.join(out_path_del, filename), sep=",", index=False)

def newton_raphson(point, nodes, GH = np.array([[-1.0, 1.0, 1.0, -1.0],[-1.0, -1.0, 1.0, 1.0]]), gh_0 = np.array([[0.0], [0.0]]), res_tol = 0.05**2, n_max = 10):
    # A newton_raphson method for the inverse mapping from Cartesian coordinates to natural
    # coordinates of a quadratic element, as used in abaqus
    # See: "Thermomechanical Modeling of Additive Manufacturing Large Parts", E Denlinger et al., Manufacturing
    # Science and Engineering, Vol 136, 2014 for details

    # point = a 3-D numpy array with vector of Cartesian coordinates for the point in question
    # nodes = 4x3 numpy array of nodal coordinates, each row is a node, ordered using the abaqus convention for quads
    # GH = natural coordinates of the nodes 1,2,3,4 following Abaqus convention for a quad
    # Note this was for ([4, 1, 5, 8]) with:
    # HR = np.array([[1.0, -1.0, -1.0, 1.0],[-1.0, -1.0, 1.0, 1.0]])
    # Will need to check if this still works
    # gh_0 = initial guess at natural coordinates for a given point, default at element centroid
    # res_tol = tolerance on the square of the residuals used to assess convergence of the Newton-Raphson
    # n_max = maximum number of iterations of the Newton-Raphson

    gh_n = gh_0 # Set at initial choice for g and h
    J_n = np.empty([3,2]) # Initialise array for storing Jacobean
    for i in range(n_max):
        # Perform the forward transformation from natural coordinates to cartesian for the current guess
        # then calculate the square-residual from the actual point Cartesian coordinates
        bases = 1.0 + gh_n*GH # 1D interpolation function in each natural coordinate
        prod_bases = np.prod(bases, axis = 0, keepdims = True)*nodes.T # 2D interpolation function
        X_i_n = np.sum(prod_bases, axis = 1)/4.0 # Cartesian coordinates for the current guess
        # calculate residual. Has to be this way round to be consistent with the minus sign in the Jacobean.
        res = point - X_i_n 
        res_L2 = np.sum(res*res) # square L2 error
        # IS stopping criterion met?
        if res_L2 < res_tol:
            i = i-1 # subtract one from i to note convergence was achieved at the previous iteration
            break
        
        # Construct the Jacobean of each Cartesian coordinate wrt g and h
        for j in range(2):
            # Minus sign consistent with residuals definition
            J_n[:,j] = -np.sum(GH[j,:]*prod_bases/bases[j,:], axis = 1)/4.0
        
        # Find the next guess of g and h using the pseduo-inverse of the non-square Jacobian
        gh_n = gh_n - linalg.pinv(J_n) @ res.reshape([3,1])

    return(gh_n, i)

def find_closest_centroid(point, centroids):
    # Find the index df element with closest centroid to point "point",
    # where "centroids" is a numpy array where each row is the centroid of an element
    # Find the distance from each point to the centroid of all elements
    cen_dist = np.sqrt(np.sum((point - centroids)*(point - centroids), axis=1))
    min_ind = np.argmin(cen_dist)
    return(min_ind)

def proj_point_on_element(point, cen, n):
    # Project point onto element with centroid "cen" and unit normal "n"
    # Calcuate vector between the point and element centroid
    v = point - cen
    # Project this vector onto the element normal to get the normal distance
    norm_dist = np.sum(v*n)
    # Move this distance along the normal such that the point lies on the surface
    point_proj = point - norm_dist*n
    return(point_proj)

def surface_mesh_from_file(file_string = [], node_file = [], el_file = []):
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
    # Create mesh_object
    mesh = SurfaceMesh(nodes, connectivity)
    
    return(mesh)

if __name__ == "__main__":
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Data_rotated", out_sub_folder="DIC_nat_coords_2")
    # node_file = '..\\outer_surface_nodes.csv'
    # el_file = '..\\outer_surface_elements.csv'
    file_string = "..\\nominal_shell_mesh_outer_surface"
    coord_labels = ["x_0_rot","y_0_rot","z_0_rot"] # list of labels for DIC coordinate labels
    # mesh = surface_mesh_from_file(node_file=node_file, el_file=el_file)
    Mesh = surface_mesh_from_file(file_string = file_string)
    # Specify that the mesh is a grid, with n_x elements in the x direction, and n_y elements in the y direction
    n_x = 84 # number of columns in the grid
    n_y = 54 # number of rows in the grid
    Mesh.define_struct_grid(n_x, n_y)
    Files.read_data()
    get_nat_coords(Files, Mesh, coord_labels=coord_labels)
    Files.dump()