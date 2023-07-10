import os
import pandas as pd
from FileSeries import *

def subtract_disp(Files):
    # Subtract the displacement from the coordinates of each dataset to retrieve the undeformed shape
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Define list of coordinate indices, and their corresponding displacement index
    column_inds = [["x", "u"],["y", "v"],["z","w"]] # index of coordinate and the corresponding displacement component
    filenames = os.listdir(Files.in_path)
    # Loop through all files in the folder
    for filename in filenames:
        src = os.path.join(Files.in_path, filename)
        # Load data from csv file
        in_data = pd.read_csv(src, sep = ",").dropna()
        # Subtract each of the displacement components from their corresponding coordinate, and 
        # create a new entry with the result
        for i, inds in enumerate(column_inds):
            # Insert entry next to the coordinates
            in_data.insert(loc = (3+i), column = inds[0]+"_0", value = (in_data[inds[0]] - in_data[inds[1]]))
    
        # Write transformed data to a new csv file
        in_data.to_csv(os.path.join(Files.out_path, filename), sep = ",", index = False)

if __name__ == "__main__":
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Data_with_outliers",out_sub_folder="Data_subtracted_disp")
    subtract_disp(Files)