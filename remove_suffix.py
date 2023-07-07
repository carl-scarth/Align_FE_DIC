# script for removing suffix (e.g. "_0.tiff") from filenames such that paraview recognises them 
# as a file series, for which the files should end in some pattern e.g. 1,2,3,4 etc

import os
import shutil
from FileSeries import *

def remove_suffix(Files, suffix="_0.tiff"):
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # suffix = suffix which is to be removed from csv files
    # Loop through each of the files in the input folder, rename, then save to output folder
    for filename in os.listdir(Files.in_path):
        src = os.path.join(Files.in_path, filename)
        dst = os.path.join(Files.out_path, filename.replace(suffix + ".csv",".csv"))
        # If writing to the same directory rename, otherwise make a copy
        if Files.in_path == Files.out_path:
            os.rename(src, dst)
        else:
            shutil.copy(src,dst)

if __name__ == "__main__":
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_1"
    folder = "..\\150kN_Data\\Camera_Pair_0_3"
    Files = FileSeries(folder=folder,in_sub_folder="100kN",out_sub_folder="")
    remove_suffix(Files)