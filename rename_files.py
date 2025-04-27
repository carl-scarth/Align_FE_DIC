# script for removing suffix (e.g. "_0.tiff") from filenames such that paraview recognises them 
# as a file series, for which the files should end in some pattern e.g. 1,2,3,4 etc

import os
import shutil
from FileSeries import *

# Directly renames files, or copies to a new folder with updated names
def rename_files(Files, suffix="_0.tiff"):
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # suffix = suffix which is to be removed from csv files
    # Loop through each of the files in the input folder, rename, then save to output folder
    [remove_suffix(File, suffix) for File in Files.files]
    for File in Files.files:
        # If writing to the same directory rename, otherwise make a copy
        if Files.in_path == Files.out_path:
            os.rename(File.src, File.dst)
        else:
            shutil.copy(File.src, File.dst)

def remove_suffix(FileObj, suffix="_0.tiff"):
    # Remove suffix (e.g. "_0.tiff") from filenames such that paraview recognises them 
    # as a file series, for which the files should end in some pattern e.g. 1,2,3,4 etc
    FileObj.out_filename = FileObj.out_filename.replace(suffix + ".csv",".csv")
    FileObj.update_dst()

if __name__ == "__main__":
    folder = "input_output"
    Files = FileSeries(folder=folder,in_sub_folder="input",out_sub_folder="output")
    rename_files(Files)