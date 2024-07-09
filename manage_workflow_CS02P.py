from FileSeries import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *

# folder = "..\\CS02P\\DIC\\Left_Camera_Pair" # Parent folder where data is located
folder = "..\\CS02P\\DIC\\Right_Camera_Pair" # Parent folder where data is located
# folder = "..\\CS04D\\DIC\\Left_Camera_Pair" # Parent folder where data is located
in_sub_folder = "Raw_Data"
out_sub_folder = "Processed_Data_Working"
# Perhaps pass instructions as a Dictionary?
rename_file = False # Do the files have a suffix which needs to be removed for paraview to recognise as a file series

# QoI = ["X","Y","Z","u","v","w","Displacement Magnitude","cu","cv","cw","exx","eyy","exy","e1","e2","Load","Crosshead","Time"]

# In general I'd rather have the loop contained outside of the functions. Think about the best way to implement that.
# Possible define a file series object with various properties? e.g. file names
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)
test_data.down_sam(8, rename_files=True)
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
# TO DO BEFORE FINISHING WITH NAT COORDS CODE:
# ADD OPTION TO DO NATURAL COORDINATE SEARCH ONLY FOR THE FIRST SAMPLE THEN DUPLICATE
# WHEN READING IN INTEGERS THESE ARE BEING CONVERTED TO FLOATS, PROBABLY WITHIN FILESERIES,
# SORT THIS OUT. MAYBE HAVE DTYPE OPTION. OR, SET DATATYPE METHOD IN FILESERIES?
# TEST METHOD FOR OUT OF BOUNDS NAT COORDS BY DELIBERATELY CREATING MESH WITH SHORT FLANGES
# CREATE ALTERNATIVE SEARCH METHOD
# TEST IF IT WORKS IN THIS ENVIRONMENT
# WRITE LIST OF INPUTS TO THE MAIN 2 FUNCTIONS. TIDY TO REMOVE MENTIONS TO EL_FOUND. BASICALLY LOOKING GOOD NOW
test_data.dump(index=True)

# Might want to print a progress bar?
