from FileSeries import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *

folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2" # Parent folder where data is located
in_sub_folder = "Raw_Data"
out_sub_folder = "Processed_Data"
# Perhaps pass instructions as a Dictionary?
rename_file = True # Do the files have a suffix which needs to be removed for paraview to recognise as a file series

QoI = ["coor.X [mm]", 
        "coor.Y [mm]", 
        "coor.Z [mm]",
        "disp.Horizontal Displacement U [mm]",
        "disp.Vertical Displacement V [mm]",
        "disp.Out-Of-Plane: W [mm]",
        "strain.Strain-global frame: Exx [ ]",
        "strain.Strain-global frame: Eyy [ ]",
        "strain.Strain-global frame: Exy [ ]",
        "strain.Maximum Principal Strain: EI [ ]",
        "strain.Minimum Principal Strain: EII [ ]"]
    
new_names = ["x", "y", "z", "u", "v", "w", "Exx", "Eyy", "Exy", "EI", "EII"]

# In general I'd rather have the loop contained outside of the functions. Think about the best way to implement that.
# Possible define a file series object with various properties? e.g. file names
test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)
test_data.read_data(sep=",")
# Generally do below, unless there is only one instruction in which case call the other function for renaming
[remove_suffix(File) for File in test_data.files]
test_data.extract_qoi(QoI, new_names=new_names, dropna = False)
subtract_disp(test_data)
# Load in a transformation matrix from file
transmat_file = os.path.join(test_data.parent_path, "image_0000_transformation_matrix.txt")
R, T = transmat_from_file(transmat_file)
# Labels of columns to which coordinate transformations are performed, and Boolean indicating which are displacements
trans_coord_labels = [["x", "y", "z"],["x_0","y_0","z_0"],["u","v","w"]]
is_displacement = [False, False, True]
transform_coords(test_data, trans_coord_labels, is_displacement = is_displacement, R=R, T=T, subscript = "rot")
test_data.dump()

# Might want to print a progress bar?
