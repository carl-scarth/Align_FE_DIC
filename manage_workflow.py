from FileSeries import *

folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_1" # Parent folder where data is located
in_sub_folder = "Raw_Data"
# Perhaps pass instructions as a Dictionary?
remove_suffix = True # Do the files have a suffix which needs to be removed for paraview to recognise as a file series

# In general I'd rather have the loop contained outside of the functions. Think about the best way to implement that.
# Possible define a file series object with various properties? e.g. file names
test_data = FileSeries(folder=folder, in_sub_folder = in_sub_folder)

# Might want to print a progress bar?
# Some NaNs to deal with tomorrow...
# Manta pair look dodgy