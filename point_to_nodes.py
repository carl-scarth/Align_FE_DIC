# Converts displacement predictions at multiple points in an element to nodes
# using least squares

import os
import pandas as pd
import numpy as np
from numpy.linalg import lstsq

def intp_coeffs(gh):
    # Returns a set of coefficients for the nodal quantities in the 
    # interpolation functions for a first-order quad elemennt given
    # natural coordinates g and h
    GH = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0]]) # Natural coordinates of the nodes
    bases = 1.0 + gh*GH
    prod_bases = np.prod(bases, axis = 1)/4.0
    return(prod_bases)

if __name__ == "__main__":
    # Directories for input/output
    folder = "E:\\MengYi_Data\\CS02P_DIC\\Left Camera Pair"
    in_subfolder = "Nat_coords_nadrop"
    out_subfolder = "Nodal_DIC"
    in_path = os.path.join(os.getcwd(),folder)
    if out_subfolder not in os.listdir(in_path):
        os.mkdir(os.path.join(in_path,out_subfolder))
    disp_inds = ["u_rot","v_rot","w_rot"]

    # Load mesh info
    node_file = "coords_undeformed.csv"
    conn_file = "element_quad.csv"
    nodes = pd.read_csv(node_file)
    conn = pd.read_csv(conn_file, dtype=int)

    # Loop over every file in the folder
    for file in os.listdir(os.path.join(in_path, in_subfolder)):
        print(file)
        src = os.path.join(in_path,in_subfolder,file)
        cloud_data = pd.read_csv(src).dropna()
        elements = sorted(cloud_data["Element"].unique().tolist())
        
        # Loop through all elements containing points and identify nodes
        # Note: Element indices are in Python convention, not Abaqus
        # Nodes are in Abaqus convention, and need to be converted to python
        node_inds = []
        [node_inds.extend((conn.iloc[element]-1).tolist()) for element in elements]
        # Sort and extract unique node ids
        node_inds = sorted(list(set(node_inds)))
        n_nodes = len(node_inds)
        
        # Sort point cloud data in terms of element number
        cloud_data = cloud_data.sort_values("Element")
        y = cloud_data[disp_inds] # Left hand side of regression, y = Xu
        
        # Extract element indices and natural coordinates for point cloud data
        el_inds = cloud_data["Element"].to_numpy()
        gh = cloud_data[["g","h"]].to_numpy()
        # Define matrix of regression coefficients, X
        X = np.zeros([cloud_data.shape[0],n_nodes])    
        for i, (el_ind, gh_point) in enumerate(zip(el_inds, gh)):
            # Get interpolation function coefficients for current point
            coeffs = intp_coeffs(gh_point)
            # Store coefficients in correct column of X, to match the node in question
            for node, coeff in zip((conn.iloc[el_ind]-1).to_numpy(), coeffs):
                X[i,node_inds == node] = coeff
        
        # Remove nodes for which there is insufficient data to reliably give displacements
        del_cnt = 0
        for i in range(len(node_inds)):
            if (X[:,i] > 0).sum()<4:
                node_inds.pop(i-del_cnt)
                del_cnt += 1
    
        X = X[:,(X > 0).sum(axis=0)>=4]
        # Also remove experimental data if there are now some rows which don't match any nodes
        y = y[(X > 0).sum(axis=1)>0]
        X = X[(X > 0).sum(axis=1)>0,:]
    
        # Solve regression y = Xu to get nodal displacements, u
        u, res, rank, s = lstsq(X, y.to_numpy(),rcond=None)
        # Store the output in the correct row for the appropriate node
        disp_nodes = np.empty((nodes.shape[0],3))
        disp_nodes.fill(np.nan)
        disp_nodes[node_inds,:] = u
        # Now save the output to csv
        disp_nodes = pd.DataFrame(disp_nodes, columns=disp_inds)
        out_frame = pd.concat((nodes, disp_nodes), axis = 1)
        out_frame.insert(0, "Node_ID",out_frame.index.values+1)
        out_frame["Load"] = cloud_data["Load"].values[0]
        print(out_frame)
        dsadsad
        dst = os.path.join(in_path,out_subfolder,file)
        out_frame.dropna().to_csv(dst,sep=",",index=False)
