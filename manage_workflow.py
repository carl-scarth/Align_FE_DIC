from FileSeries import *
from SurfaceMesh import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *
from get_nat_coords import *

# Main script for carrying out point cloud manipulations

# Define folders where input and output are to be read/written
folder = "input_output" # Parent folder
in_sub_folder = "input" # Subfolder where input point cloud data is stored
out_sub_folder = "output" # Subfolder where output is to be written
mesh_file_string = "input_output\\FE_mesh" # Identifier of node and element connectivity file

# Create list of Quantities of Interest to Extract, and define updated column names
QoI = ["x", "y", "z", "u", "v", "w", "DAQ_DAQCrossheadmm","DAQ_DAQForcekN"]
QoI_renamed = ["x", "y", "z", "u", "v", "w", "Crosshead","Force"]

# Load in transformation rotation matrix and translation vector from file
transmat_file = os.path.join(folder, "transformation_matrix.txt")
R, T = transmat_from_file(transmat_file)

# Load in mesh and create mesh object, containing nodal coordinates and connectivities
Mesh = SurfaceMesh(from_file = True, file_string=mesh_file_string)
# Specify that the mesh is an ordered grid, with n_x elements in the x direction, 
# and n_y elements in the y direction (don't do if mesh is unstructured, another
# method is used instead)
n_x = 84 # number of columns in the grid
n_y = 54 # number of rows in the grid
Mesh.define_struct_grid(n_x, n_y)

# Create FileSeries object and read data from csvs
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)
test_data.read_data(sep=",")

# Down sample across time to retain 1 in 2 of the point clouds
test_data.down_sam(2, rename_files = True)

# Extract the chosen quantities of interest, and rename columns
test_data.extract_qoi(QoI, new_names=QoI_renamed, dropna = False)

# Subtract displacements from the coordinates to back-calculate undeformed point cloud
# Append subscript "0" to column labels of output
subtract_disp(test_data, column_labels = [["x", "u"],["y", "v"],["z","w"]], subscript = "0")

# Apply the coordinate transformations specified by rotation matrix R and translation vector T
# to x,y & z, x0, y0 & z0, and u, v & w. The latter are displacements so only a rotation is applied
transform_coords(test_data, [["x", "y", "z"],["x_0","y_0","z_0"],["u","v","w"]], 
                 is_displacement = [False, False, True], R=R, T=T, subscript = "rot")

# Project DIC points onto FE model surface, and map Cartesian coordinates onto
# (isoparametric) element natural coordinates. Map "undeformed" coordinates to
# ensure compatibility with the model where the displacement may differ
# specify first_file_only = True to copy mapped coordinates from first file 
# to other increments, to save time if points are identical across increments
get_nat_coords(test_data, Mesh, in_sub = "0_rot", out_cols = [4, 9, 14, -1, -1, -1], proj_sub = "proj", first_file_only = False)

# Filter the points to retain only points between 15.0 <= z_proj <= 405.0
bound_list = ["z_proj < 15.0", "z_proj > 405.0"]
test_data.filter_by_cond(bound_list, drop_cond = True)

# Write output csvs
test_data.dump(dropna=True, index = True, index_label = "point_index")