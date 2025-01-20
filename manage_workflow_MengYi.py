from FileSeries import *
from SurfaceMesh import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *
from project_points import *
import os

# Set up FileSeries object. This contains all the code the details of
# the main directory where data is stored, sub_folder with input data,
# and subfolder to which data is written. Also code for looping over
# different files, reading and writing, and applying functions to data
folder = "E:\\MengYi_Data\\CS02P_DIC\\Right Camera Pair"
# folder = "E:\\MengYi_Data\\CS02P_DIC\\Left Camera Pair" # Parent folder of data
in_sub_folder = "Raw_Data" # Subfolder where input data is stored
out_sub_folder = "Projection_Test" # Subfolder where output data will be written (will be created if it doesn't exist already)
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)

# Use if you want to select a subset of the data (e.g. take a regular sample every 20 images)
test_data.down_sam(30, rename_files=True)

# Read in data from files 
test_data.read_data(sep=",") # Read in the data

# Use to extract specified Quantities of Interest from the data
# QoI = ["X","Y","Z","u","v","w","Displacement Magnitude","cu","cv","cw","exx","eyy","exy","e1","e2","Load","Crosshead","Time"]
#test_data.extract_qoi(QoI, new_names=new_names, dropna = False)

# Apply translations and rotations to data
# Load in a transformation matrix from file (as outputted from cloudcompare)
transmat_file = os.path.join(test_data.parent_path, "transformation_mat.txt")
R, T = transmat_from_file(transmat_file)
# Labels of columns to which coordinate transformations are performed, and Boolean indicating which are displacements
trans_coord_labels = [["x", "y", "z"],["u","v","w"]]
is_displacement = [False, True]
# Apply transformations
transform_coords(test_data, trans_coord_labels, is_displacement = is_displacement, R=R, T=T, subscript = "rot")

# Project onto mesh
# Load in mesh and create mesh object, containing nodal coordinates and
# connectivities, as well as methods for calculating centroids, normals etc
node_file = "E:\\MengYi_Data\\coords_undeformed.csv"
el_file = "E:\\MengYi_Data\\element_quad.csv"
# Construct mesh object based on connectivities, and calculate 
# element normals and centroids
FEMesh = SurfaceMesh(from_file = True, node_file=node_file, el_file=el_file)
project_points(test_data, FEMesh, in_sub = "rot")
# Write transformed data to csvs
test_data.dump()