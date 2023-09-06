# script for removing suffix (e.g. "_0.tiff") from filenames such that paraview recognises them 
# as a file series, for which the files should end in some pattern e.g. 1,2,3,4 etc

import os
import shutil
from FileSeries import *

# Directly renames files, or copies to a new folder with updated names
def remove_suffix(Files, suffix="_0.tiff"):
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # suffix = suffix which is to be removed from csv files
    # Loop through each of the files in the input folder, rename, then save to output folder
    for filename in Files.in_filenames:
        src = os.path.join(Files.in_path, filename)
        dst = os.path.join(Files.out_path, filename.replace(suffix + ".csv",".csv"))
        # If writing to the same directory rename, otherwise make a copy
        if Files.in_path == Files.out_path:
            os.rename(src, dst)
        else:
            shutil.copy(src,dst)

def remove_suffix_file_series(Files, suffix="_0.tiff"):
    # Updates the list of names in a FileSeries object to remove the suffix. To be used when multiple
    # operations are performed and a new csv is written
    Files.out_filenames = [filename.replace(suffix + ".csv",".csv") for filename in Files.in_filenames]

if __name__ == "__main__":
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Raw_data",out_sub_folder="Suffix_removed")
    remove_suffix(Files)