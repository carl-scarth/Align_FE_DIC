# Converts from x, y and z point data to natural coordinates of a finite model
# by projecting points onto the surface and then using a Newton-Raphson solver
# to determine the natural coordinates of the projected points
# Correct implementation of full search assumes the elements are arranged in a 
# regular grid. The code may return out-of-bounds natural coordinates in some cases
# if not

import pandas as pd
import numpy as np
from numpy import linalg
import warnings
from FileSeries import *
from SurfaceMesh import *
from project_points import *

def get_nat_coords(Files, Mesh, coord_labels = ["x","y","z"], in_sub = [], proj_sub = "proj", out_cols = []):
    # Read in DIC data from csv files, then apply rotations and translations
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Mesh = Mesh object containing list of nodes and elements, and methods for determining element properties
    # coord_labels = list of strings with labels of the DIC data columns containing the coordinates
    Files.apply_func_to_data(lambda x:nat_coord_search(x, Mesh, coord_labels, in_sub = in_sub, proj_sub=proj_sub),out_cols = out_cols, message = "Projecting points onto mesh and finding natural coordinates", insert_by_label = False)
    # Going to now need to add exception to FileSeries to catch wrong input combos
    # will need to update labels as kwarg rather than arg in other codes to reflect change in structure

def nat_coord_search(cloud_data, Mesh, coord_labels, in_sub = [], proj_sub = "proj"):
    # Main code for finding natural coordinates. Cloud data is the dataframe with relevent 
    # quantities for the point cloud. Mesh is a mesh object for the mesh to which the cloud
    # is being aligned.

    if in_sub:
        in_labels = ["_".join((label,in_sub)) for label in coord_labels]
    else:
        in_labels = coord_labels

    cloud_xyz = cloud_data[in_labels].to_numpy() # coordinates

    # Initialise outputs
    xyz_proj = np.empty([cloud_xyz.shape[0],3]) # xyz coordinate of projected point
    xyz_proj.fill(np.nan)
    gh = np.empty([cloud_xyz.shape[0],2]) # Array for storing the natural coordinates of each point
    el_ind = [] # list for storing element containing each point
    conv = [] # list containing how many iterations of the Newton-Raphson were required for each point
    
    # Consider moving these into another function once the inner workings of the loop have been tidied
    # Maybe something which performs the outer element search?
    for i, point in enumerate(cloud_xyz):
        # Find element with the closest centroid to the point, to use as an intial guess for the element
        min_el = find_closest_centroid(point, Mesh.centroids)
        # A for loop with a break might be more sensible... Implement on another iterations

        # Search to determine if the point is found to lie in the chosen element. If  correct element has been 
        # chosen the natural coordinates will lie between +/-1.0, if not the converged values of these coordinatess
        # indicate a direction to move in the grid until the correct element is found, or it is determined that the 
        # point does not lie over the mesh surface
        el_found = False # Boolean to check whether the element the point belongs to has been found
        gh_prev = np.array([[np.nan], [np.nan]]) # for storing previous converged value of hr when searching
        # If implementing as a for with conditional, rather than while, do this on the first iteration 
        # for an out-of-bounds element to save doing this for every point...
        while not el_found:                
            # Project point onto the surface of the current element, and store result
            point_proj = proj_point_on_element(point, Mesh.elements[min_el].centroid, Mesh.elements[min_el].n)
            xyz_proj[i,:] = point_proj
            # Determine natural coordinates for the  point using Newton-Raphson
            # Continue to update function
            (gh_i, i_conv) = newton_raphson(point_proj, Mesh.elements[min_el].nodes)
                    
            # Check if the converged natural coordinates are within the element bounds (these
            # should be +/-1.) If not, then another element is selected, and another iteration
            # of the outer while loop is performed
            if np.all((gh_i >= -1) & (gh_i <= 1)):
                # Condition for exiting the while loop
                el_found = True
            elif Mesh.is_grid:
                # Job to do on next tidy: replace the while loop with a for loop with a break.
                # Maybe set loop limit to 5??
                # There isn't a built in function which will redo the for loop iteration if necessary. 
                #  It might be tidier to do this using a while loop at the outer level, rather than a for, and
                #  controlling the incrementation. A break statement inside a conditional will prevent moving to the next iteration. Consider playing with this 
                # and see which is neater. See:
                # https://stackoverflow.com/questions/492860/python-restarting-a-loop
                # https://stackoverflow.com/questions/36573486/redo-for-loop-iteration-in-python
                # Determine position of current element in the grid to find the next element
                # Could add below to grid code?
                row_ind, col_ind = get_grid_inds(min_el,Mesh.n_x)

                # Move along to the next element in the grid based upon the converged value of the 
                # natural coordinate gives a clue of which direction in the grid the actual element lies in
                # I've corrected this to match 1,2,3,4 now. Be careful though as this will depend on in which
                # direction the h and r (g and h) increase with increasing element numbering. I think what I
                # have is a reasonable default, but it might be worth having this as an input? Perhaps orientation
                # of the grid relative to g and h is also something to think about
                if gh_i[0,0] < -1:
                    row_ind = row_ind - 1
                elif gh_i[0,0] > 1:
                    row_ind = row_ind + 1
                if gh_i[1,0] < -1:
                    col_ind = col_ind - 1
                elif gh_i[1,0] > 1:
                    col_ind = col_ind + 1
                
                # Check if search has converged, and update natural coordinates accordingly
                gh_i, el_found = update_grid_convergence(gh_i, gh_prev, row_ind, col_ind, Mesh)
                if not el_found:
                    # If neither, move to the next element in the grid
                    min_el = row_ind*Mesh.n_x + col_ind
                    gh_prev = gh_i # Store the previous natural coordinate values to track the element search
            else:
                warnings.warn("Search not yet implemented for irregular mesh, some results may be inaccurate")
                el_found = True # break out of loop
            
        # Store the index of the element in which the current point sits
        el_ind.append(min_el)
        conv.append(i_conv)
        gh[i,:] = gh_i.squeeze()
    
    # Put output in correct format
    out_frame = pd.DataFrame(xyz_proj, columns=["_".join((label, proj_sub)) for label in coord_labels])
    # Add new entries to the cloud_data dataframe
    gh = pd.DataFrame(gh, columns=["g","h"])
    out_frame = pd.concat([out_frame, gh], axis=1)
    out_frame["Element"] = pd.Series(el_ind).astype(int)
    out_frame["Conv_Iteration"] = pd.Series(conv).astype(int)
    return(out_frame)

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

def get_grid_inds(el, n_x):
    # Gets row and column index for current point in grid based on element index el,
    # and number of columns n_x
    row_ind = el//n_x # Row index. Divide by number of columns and round down to nearest integer
    col_ind = el%n_x # Remainder after the above integer division gives the column index
    return(row_ind, col_ind)

def update_grid_convergence(gh, gh_prev, row_ind, col_ind, Mesh):
    # Define grid object??? Might make this tidier
    # Check if the grid search has converged, return a logical if so and update values
    # Check if either the row or column index are outside the bounds of the grid. 
    # If so, the point doesn't belong to any element and should be deleted
    if ((row_ind < 0) | (row_ind >= Mesh.n_y) | (col_ind < 0) | (col_ind >= Mesh.n_x)):
        # Return nans for everything so the points may be deleted outside of the loop
        # min_el = np.nan - deleted to provent later issues with pandas - just get rid of probably
        # j = np.nan
        print(gh)
        print(row_ind, col_ind)
        gh = np.array([[np.nan], [np.nan]]) # Update natural coordinates
        converged = True
        print("out of bounds")

    elif np.any(((gh < -1) & (gh_prev > 1)) | ((gh > 1) & (gh_prev < -1))):
        # Check for flips from hr_n > 1 to  hr_n < -1 or vice-versa, indicating the point
        # is between the projection of the element faces onto the plane of the point, 
        # which can occur if abive a convex curve. If this happens round to the element boundary 
        if ((gh[0,0] < -1) | (gh[0,0] > 1)):
            gh[0,0] = round(gh[0,0])
        if ((gh[1,0] < -1) | (gh[1,0] > 1)):
            gh[1,0] = round(gh[1,0])
        converged = True # Break out of the while loop
    else:
        # If neither, need to continue the search 
        converged = False

    return(gh, converged)

if __name__ == "__main__":
    # Create file series and load in data
    folder = "E:\\MengYi_Data\\CS02P_DIC\\Right Camera Pair"
    Files = FileSeries(folder=folder,in_sub_folder="Trimmed_Data", out_sub_folder="Nat_coords_Test")

    # Load in mesh and create mesh object, containing nodal coordinates
    # and connectivities, as well as methods for calculating centroids,
    # normals etc
    node_file = "E:\\MengYi_Data\\coords_undeformed.csv"
    el_file = "E:\\MengYi_Data\\element_quad.csv"
    # Construct mesh object based on connectivities, and calculate 
    # element normals and centroids
    Mesh = SurfaceMesh(from_file = True, node_file=node_file, el_file=el_file)

    # Use if mesh is ordered as a structured grid - as this info can help improve
    # the results if a point is initially assigned to the wrong element
    # Specify that the mesh is a grid, with n_x elements in the x direction, and n_y elements in the y direction
    #n_x = 84 # number of columns in the grid
    #n_y = 54 # number of rows in the grid
    #Mesh.define_struct_grid(n_x, n_y)

    # Project DIC onto mesh surface and determine natural coordinates
    Files.read_data()
    out_cols = [3, 6, 9, -1, -1, -1, -1] # Index of colums where output data is to be inserted (will need to subtract one when using in main file - pandas has read in the index here...)
    get_nat_coords(Files, Mesh, in_sub = "rot", out_cols=out_cols)
    Files.dump()
