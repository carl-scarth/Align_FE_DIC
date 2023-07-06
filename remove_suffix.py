# script for removing suffix (e.g. "_0.tiff") from filenames such that paraview recognises them 
# as a file series, for which the files should end in some pattern e.g. 1,2,3,4 etc

import os

def remove_suffix(dataset="Raw_Data", suffix="_0.tiff", out_folder = ""):
    # dataset = subfolder in which csv files are contained
    # suffix = suffix which is to be removed from csv files
    # Path to parent folder of data
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    wd = os.path.abspath("") # Get working directory
    in_path = os.path.join(wd,folder,dataset)
    # Possibly pass in kwargs?
    if out_folder == "":
        out_path = in_path
    
    # Try to take this outside of function
    for filename in os.listdir(in_path):
        src = os.path.join(in_path, filename)
        dst = os.path.join(out_path, filename.replace(suffix + ".csv",".csv"))
        os.rename(src, dst)

if __name__ == "__main__":
    remove_suffix()