
# Loads in DIC data and subtracts the displacement from each set of coordinates to retrieve 
# the undeformed geometry.

import os
import pandas as pd
from FileSeries import *

def subtract_disp(Files):
    # Subtracts the displacement from the coordinates of each dataset to retrieve the undeformed shape
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Define list of coordinate indices, and their corresponding displacement index
    column_inds = [["x", "u"],["y", "v"],["z","w"]] 
    filenames = os.listdir(Files.in_path)
    # Loop through all files in the folder
    for filename in filenames:
        print(filename)
        src = os.path.join(Files.in_path, filename)
        # Load data from csv file
        in_data = pd.read_csv(src, sep = ",").dropna()
        # Subtract each of the displacement components from their corresponding coordinate, and 
        # create a new entry with the result
        for inds in column_inds:
            in_data[inds[0]+"_0"] = in_data[inds[0]] - in_data[inds[1]]

        # Try to use the below code to insert the columns at a specidied index, rather than appending and shuffling
        in_data.insert()
        new_cols = {"x_proj": [0,6], "y_proj": [1,7],"z_proj":[2,8]} 
        [cloud_data.insert(loc=value[1], column = key, value = pd.Series(xyz_proj[:,value[0]])) for key, value in new_cols.items()]    
    # Re-shuffle the order of the columns to place the undeformed coordinates next to the deformed ones
    new_order = ["x", "y", "z", "x_0", "y_0", "z_0", "u", "v", "w", "corrU", "corrV", "corrW", "e1", "e2"]
    out_data = in_data[new_order]
    
    # Write transformed data to a new csv file
    out_data.to_csv(os.path.join(out_path, filename), sep = ",", index = False)

if __name__ == "__main__":
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Data_with_outliers",out_sub_folder="Data_subtracted_disp")
