from FileSeries import *
from SurfaceMesh import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *
from get_nat_coords import *
from trim_bounds import *

# folder = "..\\CS02P\\DIC\\Left_Camera_Pair" # Parent folder where data is located
folder = "..\\CS02P\\DIC\\Right_Camera_Pair" # Parent folder where data is located
# folder = "..\\CS04D\\DIC\\Left_Camera_Pair" # Parent folder where data is located
in_sub_folder = "Raw_Data_Subset"
out_sub_folder = "Processed_Data_Working"
# Perhaps pass instructions as a Dictionary?
rename_file = False # Do the files have a suffix which needs to be removed for paraview to recognise as a file series

# QoI = ["X","Y","Z","u","v","w","Displacement Magnitude","cu","cv","cw","exx","eyy","exy","e1","e2","Load","Crosshead","Time"]

# In general I'd rather have the loop contained outside of the functions. Think about the best way to implement that.
# Possible define a file series object with various properties? e.g. file names
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)
# test_data.down_sam(8, rename_files=True)
test_data.read_data(sep=",")
# Generally do below, unless there is only one instruction in which case call the other function for renaming
#test_data.extract_qoi(QoI, new_names=new_names, dropna = False)
# subtract_disp(test_data, column_labels = [["x", "u"],["y", "v"],["z","w"]])

# Load in a transformation matrix from file
transmat_file = os.path.join(test_data.parent_path, "rotation_matrix.txt")
R, T = transmat_from_file(transmat_file)
# Labels of columns to which coordinate transformations are performed, and Boolean indicating which are displacements
trans_coord_labels = [["x", "y", "z"],["u","v","w"]]
is_displacement = [False, True]
transform_coords(test_data, trans_coord_labels, is_displacement = is_displacement, R=R, T=T, subscript = "rot")

# Trim data which is clearly out-of-bounds
bound_list = ["z_rot < 8.4 & y_rot < -2.15 & x_rot > 59.8",
              "x_rot > 61.4 & y_rot < -3.5 & z_rot > 410.0"]
test_data.filter_by_cond(bound_list, drop_cond = True, dropna = False, na_col = "Time")

# Load in mesh and create mesh object, containing nodal coordinates
# and connectivities, as well as methods for calculating centroids,
# normals etc
mesh_file_str = "..\\new_spar_mesh_outer_surface"
# Construct mesh object based on connectivities, and calculate 
# element normals and centroids
mesh = SurfaceMesh(from_file = True, file_string=mesh_file_str)

# Specify that the mesh is a grid, with n_x elements in the x direction, and n_y elements in the y direction
n_x = 84 # number of columns in the grid
n_y = 54 # number of rows in the grid
mesh.define_struct_grid(n_x, n_y)

# Define index of colums where output data is to be inserted. (Can leave empty to append to end, or specify each column with -1)
out_cols = [2, 5, 8, -1, -1, -1]
get_nat_coords(test_data, mesh, in_sub = "rot", out_cols=out_cols, first_file_only = True)

# Trim data which is outside of the required window
# Define list of bounds, outside of which to trim data
bound_list = ["y_proj > 145", # points in the far radius of the ends,
              "z_proj >= 185.0 & z_proj <= 235.0 & y_proj > 138.75", # Points in the far radius in the central region
              "z_proj >= 60.0 & z_proj < 185.0 & y_proj > -z_proj / 20.0 + 148.0", # Point in the far radius in the ramp
              "z_proj > 235.0 & z_proj <= 360.0 & y_proj > z_proj / 20.0 + 127.0"]

# trim_bounds(test_data, bound_list, drop_cond=True, dropna=False, na_col = "Time")
test_data.filter_by_cond(bound_list, drop_cond = True, dropna = False, na_col = "Time")

# Finally test on Meng Yi's dataset and email
test_data.dump(index=True, dropna = True)
test_data.dump_nas()

# Might want to print a progress bar?
# ALSO ADD OPTION OF PROCESSING FILES 1 AT A TIME TO SAME MEMORY