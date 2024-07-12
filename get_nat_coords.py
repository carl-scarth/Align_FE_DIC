# Converts from x, y and z point data to natural coordinates of a finite model
# by projecting points onto the surface and then using a Newton-Raphson solver
# to determine the natural coordinates of the projected points

import pandas as pd
import numpy as np
from numpy import linalg
from FileSeries import *
from SurfaceMesh import *
from project_points import *
from interpolate_data import intp_nodes_to_cloud

def get_nat_coords(Files, Mesh, coord_labels = ["x","y","z"], in_sub = [], proj_sub = "proj", out_cols = [], first_file_only = False):
    # Read in DIC data from csv files, then apply rotations and translations
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Mesh = Mesh object containing list of nodes and elements, and methods for determining element properties
    # coord_labels = list of strings with labels of the DIC data columns containing the coordinates
    # in_sub = subscript on coordinate labels used to identify input columns
    # proj_sub = subscript applied to projected coordinate output column labels 
    # out_cols = list of integer indices for natural coordinate output. IF not specified added to end of dataframe
    Files.apply_func_to_data(lambda x:nat_coord_search(x, Mesh, coord_labels, in_sub = in_sub, proj_sub=proj_sub),out_cols = out_cols, message = "Projecting points onto mesh and finding natural coordinates", insert_by_label = False, first_file_only=first_file_only)

def nat_coord_search(cloud_data, Mesh, coord_labels, in_sub = [], proj_sub = "proj", n_iter = 4):
    # Main code for finding natural coordinates. 
    # Cloud data = dataframe containing point cloud data
    # Mesh = mesh object onto which the cloud is to be mapped
    # coord_labels = string used to identify datframe columns containing input coordinates
    # proj_sub = subscript added to coordinates of projected points in output dataframe
    # n_iter = maximum number of search iterations. Default 4 is worst case for if point is on a corner

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
    # Loop over all points and perform the natural coordinate search for each point
    for i, point in enumerate(cloud_xyz):
        # Find element with the closest centroid to the point, to use as an intial guess for the element.
        # Returns the first n_iter elements
        nearest_els = find_closest_centroid(point, Mesh.centroids, sorted_points = True, n_sort=n_iter)
        min_el = nearest_els[0]
        for j in range(n_iter):
            # Project point onto the surface of the current element, and store result
            point_proj = proj_point_on_element(point, Mesh.elements[min_el].centroid, Mesh.elements[min_el].n)
            xyz_proj[i,:] = point_proj
            # Determine natural coordinates for the  point using Newton-Raphson
            gh_i, *_ = newton_raphson(point_proj, Mesh.elements[min_el].nodes)

            # Check if the converged natural coordinates are within the element bounds
            # (should be +/-1.) If not, move to the next element in the search.
            if np.all((gh_i >= -1.0) & (gh_i <= 1.0)):
                # Exit the loop
                break
            elif Mesh.is_grid:
                # Use the grid search method if mesh is a structured grid
                if j == 0:
                    gh_prev = np.array([[np.nan], [np.nan]]) # for tracking previous converged value of gh

                # Find the next element to consider in the search, and check if stopping criteria are met
                gh_i, gh_prev, min_el, stop, update_coords = update_grid_search(gh_i, gh_prev, min_el, Mesh)
                # Break out of loop if stopping criteria met
                if stop:
                    # Update projected coordinates to reflect updated element and natural coordinates
                    if update_coords:
                        xyz_proj[i,:] = intp_nodes_to_cloud([0], gh_i.T, Mesh.nodes, Mesh.elements[min_el].connectivity.reshape((1,-1)), GH = [], skip_nodes = 0).squeeze()
                    break
                elif j==n_iter-1:
                    # Catch to remove non-feasible points if no other stopping criteria are met in the final loop
                    gh_i = np.array([[np.nan], [np.nan]])

            else:
                # Otherwise update search based upon element with next closest centroid
                if j == 0:
                    res_all = []
                    gh_all = []
                    min_el_all = []
                gh_i, gh_all, res_all, min_el, min_el_all = update_cen_search(j, n_iter, nearest_els, gh_i, gh_all, min_el, min_el_all, res_all)
                if j == n_iter-1:
                    #  Update projected coordinates to reflect updated element and natural coordinates
                    xyz_proj[i,:] = intp_nodes_to_cloud([0], gh_i.T, Mesh.nodes, Mesh.elements[min_el].connectivity.reshape((1,-1)), GH = [], skip_nodes = 0).squeeze()
            
        # Store the index of the element in which the current point sits
        el_ind.append(min_el)
        gh[i,:] = gh_i.squeeze()
    
    # Put output in correct format
    out_frame = pd.DataFrame(xyz_proj, columns=["_".join((label, proj_sub)) for label in coord_labels])
    # Add new entries to the cloud_data dataframe
    gh = pd.DataFrame(gh, columns=["g","h"])
    out_frame = pd.concat([out_frame, gh], axis=1)
    out_frame["Element"] = pd.Series(el_ind).astype(int)
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

def update_grid_search(gh, gh_prev, el_ind, Mesh, update_coords = False):
    # If the natural coordinates have converged to non-feasible values, find the next element
    # in the search, and check if stopping criteria have been met 
    # Stopping criteria are either:
    # i) search has moved out of bounds of the mesh
    # ii) the point lies between two elements
    
    # Find position of current element in the grid
    row_ind, col_ind = get_grid_inds(el_ind, Mesh.n_x)
    # Move to the next element in the grid based upon the converged natural coordinate
    # values, g > 1 means move right in grid, h > 1 move up and vice versa for values < -1.  
    if gh[0,0] < -1:
        row_ind = row_ind - 1
    elif gh[0,0] > 1:
        row_ind = row_ind + 1
    if gh[1,0] < -1:
        col_ind = col_ind - 1
    elif gh[1,0] > 1:
        col_ind = col_ind + 1
                
    # Check if stopping criteria met, and update natural coordinates accordingly
    if ((row_ind < 0) | (row_ind >= Mesh.n_y) | (col_ind < 0) | (col_ind >= Mesh.n_x)):
        # Return nans for everything so the points may be deleted outside of the loop
        print(gh)
        print(row_ind)
        print(col_ind)
        gh = np.array([[np.nan], [np.nan]]) # Update natural coordinates
        stop = True
        print("out of bounds")
        asdsadsad # I haven't yet had a mesh to check this works (try to find one???)

    elif np.any(((gh < -1) & (gh_prev > 1)) | ((gh > 1) & (gh_prev < -1))) and (np.all(np.abs(gh)<1.2)):
        # Check for flips from gh_n > 1 to  gh_n < -1 or vice-versa, indicating the point
        # is between two elements which can occur if above a convex curve. 
        # If this happens round to either element boundary
        if ((gh[0,0] < -1) | (gh[0,0] > 1)):
            gh[0,0] = round(gh[0,0])
        if ((gh[1,0] < -1) | (gh[1,0] > 1)):
            gh[1,0] = round(gh[1,0])
        update_coords = True
        stop = True
    else:
        # If neither, move to the next element in the grid
        el_ind = row_ind*Mesh.n_x + col_ind
        gh_prev = gh # Store the previous natural coordinate values to track the element search
        stop = False

    return(gh, gh_prev, el_ind, stop, update_coords)

def update_cen_search(j, n_iter, nearest_els, gh, gh_all, el_ind, el_ind_all, res_all, out_bounds_tol = 0.05):
    # If the natural coordinates have converged to non-feasible values, find the next element
    # in the search, based upon element with next closest centroid. Used if the mesh is not a 
    # regular grid. Stopping criterion is dictated external to function, when for loop ends
                
    # Store sample residual. Note the residual from Newton-Raphson isn't reliable as the 
    # projected position moves with different element guesses, so use the maximum extent to 
    # which each natural coordinate is outside the feasible domain
    res_all.append(max(np.sum(-gh[gh < -1]-1.0), np.sum(gh[gh > 1]-1.0)))
    el_ind_all.append(el_ind) # Store current element guess
    gh_all.append(gh) # Store gh for current guess

    # Move to the next element in the search
    if j < n_iter-1:
        el_ind = nearest_els[j+1]   # Update element number
    else:
        # On final iteration choose element resulting in lowest residual. Detect if this 
        # residual is too high and point is outside mesh. If not, and result is still outside 
        # of feasible values, round to 1.0 or -1.0 to give point on element boundary.
        min_ind = [ind for ind, min_res in enumerate(res_all) if min_res == min(res_all)][0]
        gh = gh_all[min_ind]
        el_ind = el_ind_all[min_ind]
        # Detect if element is outside of mesh using maximum tolerance infeasible natural 
        # coordinate values. Default value more than 1/20th of element.
        if np.any((gh-1.0>out_bounds_tol) | (gh+1.0<-out_bounds_tol)):
            # Return nans for everything so the points may be deleted outside of the loop
            gh = np.array([[np.nan], [np.nan]])
        # Round results if necessary
        gh[gh>1.0] = 1.0
        gh[gh<-1.0] = -1.0
    
    return(gh, gh_all, res_all, el_ind, el_ind_all)


if __name__ == "__main__":
    # Create file series and load in data
    folder = "..\\CS02P\\DIC\\Right_Camera_Pair"
    Files = FileSeries(folder=folder,in_sub_folder="Processed_Data_Working", out_sub_folder="Nat_coords")

    # Load in mesh and create mesh object, containing nodal coordinates
    # and connectivities, as well as methods for calculating centroids,
    # normals etc
    file_string = "..\\new_spar_mesh_outer_surface"
    # Construct mesh object based on connectivities, and calculate 
    # element normals and centroids
    Mesh = SurfaceMesh(from_file = True, file_string=file_string)

    # Specify that the mesh is a grid, with n_x elements in the x direction, and n_y elements in the y direction
    n_x = 84 # number of columns in the grid
    n_y = 54 # number of rows in the grid
    Mesh.define_struct_grid(n_x, n_y)

    # Project DIC onto mesh surface and determine natural coordinates
    Files.read_data(dropna=True)
    # Covert index to integer
    Files.update_datatype('int',index=[0])
    out_cols = [3, 6, 9, -1, -1, -1] # Index of colums where output data is to be inserted (will need to subtract one when using in main file - pandas has read in the index here...)
    get_nat_coords(Files, Mesh, in_sub = "rot", out_cols=out_cols, first_file_only = True)
    Files.dump()
