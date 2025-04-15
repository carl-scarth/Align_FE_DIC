 
# Takes in test data, extracts quantities of interest, then dumps resulting csvs into a folder

from FileSeries import *

def extract_qoi(Files, QoI, new_names = []):
    # Extracts columns of interest from a csv files, then removes rows with missing values
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # QoI = list of strings indicating which columns are to be extracted
    # new_names = list of updated names which are to be applied to the columns. Must be the same length as QoI
    
    # Loop through all files in the folder
    Files.read_data(sep=",")
    Files.extract_qoi(QoI, new_names)
    # Write selected data to paths
    Files.dump()

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

    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Suffix_removed",out_sub_folder="Extracted_qoi")
    extract_qoi(Files, QoI, new_names=new_names)
    Files.dump()