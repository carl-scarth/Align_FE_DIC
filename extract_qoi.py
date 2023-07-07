
# Takes in test data, extracts quantities of interest, then dumps resulting csvs into a folder

import os
import pandas as pd
import warnings

def extract_qoi(QoI, folder = "", in_sub_folder = "Raw_Data", out_sub_folder = "Data_with_outliers", new_names = []):
    # Extracts columns of interest from a csv files, then removes rows with missing values
    # dataset = subfolder in which csv files are contained
    # suffix = suffix which is to be removed from csv files
    # Path to parent folder of data
    
    wd = os.path.abspath("") # Get working directory
    in_path = os.path.join(wd,folder,in_sub_folder)
    # Possibly pass in kwargs?
    if out_sub_folder == "":
        out_path = in_path
    else:
        out_path = os.path.join(wd,folder,out_sub_folder)
        
    # Create out_path if it doesn't already exist
    if not(os.path.isdir(out_path)):
        os.mkdir(out_path)
        
    # Get filenames
    filenames = os.listdir(in_path)
    # Loop through all files in the folder
    for filename in filenames:
        print(filename)
        src = os.path.join(in_path, filename)
        # Load data from csv file, filter for QoI, then drop NaN values
        in_data = pd.read_csv(src, sep = ",")
        out_data = in_data.filter(items=QoI).dropna()
        # Rename Quantities of Interest if required
        if new_names:
            if len(out_data.columns) == len(new_names):
                out_data.columns = new_names
            else:
                warnings.warn("Specified number of new_names did not match the selected number of Colummns. Could not rename.")

        # Write selected data to paths
        out_data.to_csv(os.path.join(out_path, filename), sep = ",", index = False)

if __name__ == "__main__":
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

    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    extract_qoi(QoI, folder = folder, new_names=new_names)