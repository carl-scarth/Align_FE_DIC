from FileSeries import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *

# Set up FileSeries object. This contains all the code the details of
# the main directory where data is stored, sub_folder with input data,
# and subfolder to which data is written. Also code for looping over
# different files, reading and writing, and applying functions to data
# folder = "E:\\MengYi_Data\\CS02P_DIC\\Right Camera Pair"
folder = "E:\\MengYi_Data\\CS02P_DIC\\Left Camera Pair" # Parent folder of data
in_sub_folder = "Raw_Data" # Subfolder where input data is stored
out_sub_folder = "Transformed_Data" # Subfolder where output data will be written (will be creating if it doesn't exist already)
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)

# Use if you want to select a subset of the data (e.g. take a regular sample every 20 images)
test_data.down_sam(20, rename_files=True)

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

# Write transformed data to csvs
test_data.dump()