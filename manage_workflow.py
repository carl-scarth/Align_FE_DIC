from FileSeries import *
from rename_files import *

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
# Add filter function here
test_data.dump()

# Might want to print a progress bar?
