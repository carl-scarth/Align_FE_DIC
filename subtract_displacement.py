from FileSeries import *
import numpy as np

def subtract_disp(Files, column_labels = [["x", "u"],["y", "v"],["z","w"]]):
    # Subtract the displacement from the coordinates of each dataset to retrieve the undeformed shape
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Define list of coordinate indices, and their corresponding displacement index
    # Column inds = list of coordinate labels and the corresponding displacement component
    Files.apply_func_to_data(lambda x:subtract_columns(x,column_labels), [labels[0] for labels in column_labels], "0", message = "Subtracting displacement")

def subtract_columns(data, column_labels = [["x", "u"],["y", "v"],["z","w"]]):

    new_cols = np.empty((data.shape[0],len(column_labels)))    
    # for labels in column_labels:
    for i, labels in enumerate(column_labels):
        new_cols[:,i] = (data[labels[0]] - data[labels[1]]).to_numpy()

    return(new_cols)

if __name__ == "__main__":
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Extracted_qoi",out_sub_folder="Data_subtracted_disp_2")
    Files.read_data()
    subtract_disp(Files)
    Files.dump()
