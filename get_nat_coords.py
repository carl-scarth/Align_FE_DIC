# Converts from x, y and z point data to natural coordinates of a finite model
# by projecting points onto the surface and then using a Newton-Raphson solver
# to determine the natural coordinates of the projected points
# Assumes the finite elements are arranged in a regular grid

# This will be better packaged into different functions etc
# Suggest that a function for element operations, and those
# only on the nodes would be best

# Certainly package some of the Newton Raphson element search functions to 
# Make the different loop/conditional statements less confusing to read

import numpy as np
from numpy import linalg
import pandas as pd
import os

# Names of all of the input files and folders
in_path = os.path.abspath("") + "\\Data_rotated\\"
out_path = os.path.abspath("") + "\\Data_nat_coords\\" # Path used to store processed outputs
out_path_del = os.path.abspath("") + "\\Data_proj_del\\" # Path used to store processed outputs
# Check if out_path and out_path_del exist, and if not create them
if not(os.path.isdir(out_path)):
    os.mkdir(out_path)
if not(os.path.isdir(out_path_del)):
    os.mkdir(out_path_del)

node_file = 'outer_surface_nodes.csv'
element_file = 'outer_surface_elements.csv'

# Read element connectivities (stick with numpy as better for linear algebra)
connectivity = np.loadtxt(element_file, dtype = int, delimiter = ',', skiprows = 1)
n_el = connectivity.shape[0]
# Re-write node indices from abaqus to python indexing
connectivity = connectivity - 1

# read nodal coordinates
nodes = np.loadtxt(node_file, delimiter = ',', skiprows = 1)

# Define grid dimensions - in the long run it would be nice to read this in from
# the model somehow if possible
n_x = 84 # number of columns in the grid
n_y = 54 # number of rows in the grid

# To project onto the surface we now need to calculate element normals.
# First extract the nodes for each element.
surf_nodes = np.empty([n_el,3,4])
for i, element in enumerate(connectivity):
    surf_nodes[i,:,:] = nodes[element,:].T

# Determine two unit vectors on the element surface
a = surf_nodes[:,:,1] - surf_nodes[:,:,0]
b = surf_nodes[:,:,3] - surf_nodes[:,:,0]
a = a/np.sqrt(np.sum(a*a, axis=1)).reshape([n_el,1])
b = b/np.sqrt(np.sum(b*b, axis=1)).reshape([n_el,1])
# Get unit normal by taking the cross product of a and b
n = np.cross(a,b)

# Calculate the element centroids
cen = np.sum(surf_nodes,2)/4.0

# Define natural coordinates of the selected nodes ([4, 1, 5, 8]). Restrict this to one face of the the brick
HR = np.array([[1.0, -1.0, -1.0, 1.0],[-1.0, -1.0, 1.0, 1.0]])
# Define initial guess for h and r
hr_0 = np.array([[0.0], [0.0]])
# tolerance on the residuals used to assess convergence of the Newton-Raphson
# (given as the square of the distance to avoid doing the sqrt for each iteration)
res_tol = 0.05**2 
n_max = 10 # maximum number of iterations of the Newton-Raphson

filenames = os.listdir(in_path)
# Loop through all point cloud files in the folder
for filename in filenames:
    print(filename)

    # Read in point cloud
    src = os.path.join(in_path, filename)
    cloud_data = pd.read_csv(src, sep = ",").dropna()
    cloud_xyz = cloud_data[["x_0","y_0","z_0"]].to_numpy() # coordinates
    # numpy code if still needed
    # cloud_data = np.loadtxt(cloud_file,delimiter= ',', skiprows = 1)
    # Get csv headers
    #with open(cloud_file,'r') as f:
        #cloud_header = f.readline().strip()

    # cloud_xyz = cloud_data[:,0:3] # coordinates
    n_cloud = cloud_xyz.shape[0]
    # Here we now need to find the distance of all points in the cloud to the element centroid.
    
    # Loop through each of the points in the cloud and project onto the surface of the mesh. 
    # Perform Newton-Raphson to then determine natural coordinates of this point
    xyz_proj = np.empty([n_cloud,3]) # xyz coordinate of projected point
    hr = np.empty([n_cloud,2]) # Array for storing the natural coordinates of each point
    el_ind = [] # list for storing element containing each point
    conv = [] # list containing how many iterations of the Newton-Raphson were required for each point
    for i, point in enumerate(cloud_xyz):
        # Find the distance from each point to the centroid of all elements
        cen_dist = np.sqrt(np.sum((point - cen)*(point - cen), axis=1))
    
        # As an initial guess, assume that the point lies in the element with
        # the closest centroid
        min_el = np.argmin(cen_dist)
        # Boolean to check whether the element the point belongs to has been found
        el_found = False
        hr_n_old = np.array([[np.nan], [np.nan]]) # for storing previous converged value of hr when searching
        # If implementing as a for with conditional, rather than while, do this on the first iteration 
        # for an out-of-bounds element to save doing this for every point...
    
        # Project the point onto the element surface and find the natural coordinates of the projected
        # point in the coordinate system of that element. Convergence to coordinates outside of the 
        # bounds of the element indicates that the incorrect element has been chosen. These
        # then dictate a direction to move in the grid, and a new element is chosen until the correct
        # element is found, or it is determined that the point does not belong to any element

        # A for loop with a break might be more sensible...
        while not el_found:
            # Calcuate Vector between centroid of closest element and the current point
            v = point - cen[min_el,:]

            # Project this vector onto the normal of the closest element to get the normal distance
            norm_dist = np.sum(v * n[min_el,:])
            # Move this distance along the normal such that the point lies on the surface
            point_proj = point - norm_dist*n[min_el,:]
            xyz_proj[i,:] = point_proj
        
            # Determine natural coordinates for the  point using Newton-Raphson
            hr_n = hr_0 # Set initial choice of h and r
            # (I don't think the below line is needed - dates back to when this was a while loop
            # delete if the code works without it)
            # res_L2 = 1e8 # Reset value of the square of the residuals
            for j in range(n_max):
                bases = 1.0 + hr_n*HR
                prod_bases = np.prod(bases, axis = 0) * surf_nodes[min_el,:,:]
                X_i_n = np.sum(prod_bases, axis = 1)/4.0
                # calculate residual
                # Note that I'm trying to find the root of the residual, = point_proj - X_i_n 
                # It has to be this way round to be consistent with the minus sign in calculation of the Jacobean.
                # If flipping the residual definition the other way round the this minus sign needs to be removed also
                res = point_proj - X_i_n 
                # Work out the square of the euclidean norm of the residual
                res_L2 = np.sum(res*res)
                # If the stopping criterion is met the loop may be terminated
                if res_L2 < res_tol:
                    j = j-1 # subtract one from j to note iteration on which convergence was achieved
                    break
            
                # Construct the Jacobean of each coordinate wrt h and r
                J_n = np.empty([3,2]) # Consider initialising this outside of the loop
                for k in range(2):
                    J_n[:,k] = -np.sum(HR[k,:]*prod_bases/bases[k,:], axis = 1)/4.0

                # Perform the next step of the Newton-Raphson using the pseduo-inverse
                # of the Jacobian, as this is non-square
                hr_n = hr_n - linalg.pinv(J_n) @ res.reshape([3,1])

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
                if hr_n[0,0] < -1:
                    row_ind = row_ind + 1
                elif hr_n[0,0] > 1:
                    row_ind = row_ind - 1
                if hr_n[1,0] < -1:
                    col_ind = col_ind - 1
                elif hr_n[1,0] > 1:
                    col_ind = col_ind + 1
            
                # Check to see if either the row or column index are outside the bounds of the grid. 
                # If so, the point doesn't belong to any element and should be deleted
                if ((row_ind < 0) | (row_ind >= n_y) | (col_ind < 0) | (col_ind >= n_x)):
                    # Return nans for everything so the points may be deleted outside of the loop
                    min_el = np.nan
                    j = np.nan
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

    # Remove points which were found to lie outside of the mesh, identified using nans
    nan_ind = np.isnan(hr[:,0])
    hr = hr[np.logical_not(nan_ind),:]
    xyz_del = xyz_proj[nan_ind,:]
    xyz_proj = xyz_proj[np.logical_not(nan_ind),:]
    cloud_data_del = cloud_data[nan_ind]
    cloud_data = cloud_data[~nan_ind]
    el_ind = [el for el in el_ind if not np.isnan(el)]
    conv = [n for n in conv if not np.isnan(n)] 
    # write adjusted point cloud, and deleted points to csv files
    # Dictionary of new column names and their positions in the array and output dataframe respectively
    new_cols = {"x_proj": [0,6], "y_proj": [1,7],"z_proj":[2,8]} 
    [cloud_data.insert(loc=value[1], column = key, value = pd.Series(xyz_proj[:,value[0]])) for key, value in new_cols.items()]
    [cloud_data_del.insert(loc=value[1], column = key, value = pd.Series(xyz_del[:,value[0]])) for key, value in new_cols.items()]
    # Also append converged, natural coordinates, element index and  number of iterations
    cloud_data = pd.concat([cloud_data, pd.DataFrame(hr, columns=["h", "r"])],axis=1)
    cloud_data["Element"] = el_ind
    cloud_data["Convergence Iteration"] = conv
    cloud_data.to_csv(os.path.join(out_path, filename), sep=",", index=False)
    cloud_data_del.to_csv(os.path.join(out_path_del, filename), sep=",", index=False)

