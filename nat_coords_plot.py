# Load in a list of element numbers and natural coordinates then convert for plotting on a regular grid

import numpy as np
import pandas as pd
import os

folder = "..\\CS02P\\DIC\\Left_Camera_Pair\\Interpolated_Data_200kN"
file_str = 'Image_Inc_16'

# Load data from csvs
# data = np.loadtxt(os.path.join(folder,file_str+'.csv'),skiprows=1,delimiter=',')
data = pd.read_csv(os.path.join(folder,file_str+'.csv'))
elements = data["Element"].to_numpy().astype(int)
gh = data[["g","h"]].to_numpy()
data = data[["u_rot","v_rot","w_rot"]]

# Define grid dimensions
n_x = 84 # number of columns in the grid
n_y = 54 # number of rows in the grid

# Get row and column index for each point using the element number

row_ind = elements//n_x # Row index. Divide by number of columns and round down to nearest integer
col_ind = elements%n_x # Remainder after the above integer division gives the column index

# This will change depending upon the orientation of the model. I've followed Jean's convention
# This will be easier if I have a more sensible definintion of h...
g_plot = (row_ind*2 + 1 - gh[:,0]) 
h_plot = (col_ind*2 + gh[:,1] + 1)

out_frame = pd.DataFrame({"x_plot":h_plot, "y_plot":g_plot, "z_plot" : np.zeros(g_plot.shape)})
out_frame = pd.concat((out_frame, data),axis=1)
out_frame.to_csv("nat_coords_plot.csv",sep=",",index=False)