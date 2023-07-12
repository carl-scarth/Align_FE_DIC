import os
import pandas as pd
import numpy as np
from FileSeries import *

def transform_coords(Files, R=[], T=[]):
    # Read in DIC data from csv files, then apply rotations and translations as specified
    # via a text file
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    
    # Loop through all files in the input folder and perform the transformation
    for filename in Files.in_filenames:
        print(filename)
        src = os.path.join(Files.in_path, filename)
        # Load data from csv file
        in_data = pd.read_csv(src, sep = ",").dropna()

        # Perform the new rotation for different sets of inputs, if they exist, and insert transformed data
        # next to original in the dataframe
        # Consider passing with arguments via a dictionary when sorting out the high level workflow
        insert_columns_with_func(in_data, ["x","y","z"], "_rot", lambda xyz:rotate_translate(xyz,R,T))
        insert_columns_with_func(in_data, ["x_0","y_0","z_0"], "_rot", lambda xyz:rotate_translate(xyz,R,T))
        insert_columns_with_func(in_data, ["u","v","w"], "_rot", lambda xyz:rotate_translate(xyz,R))        
    
        # Write transformed data to a new csv file
        in_data.to_csv(os.path.join(Files.out_path, filename), sep = ",", index = False)

def rotate_translate(xyz, R=[], T=[]):
    # apply a coordinate rotation and transformation to array of coordinates, where xyz is a N x d (N = number of 
    # samples, d is dimension of coordinates) array of coordinates, R is a d x d rotation matrix, and T is a 
    # d-dimensional translation vector. If either R or T are not passed the corresponding transformation is not performed
    if len(R)!=0:
        xyz = (R@xyz.T).T
    if len(T)!=0:
        xyz = xyz + T

    return(xyz)

def transmat_from_file(trans_mat_file, delimiter=" ", header = False):
    # Load in transformation matrix from a text file
    # trans_mat_file = file where transformation matrix is defined. This is in the format outputted by CloudCompare
    # where the first 3 columns x 3 rows is a rotation matrix, and the 4th column with 3 rows is the translation vector
    # Standard format is delimited by spaces with no header, but can pass alternatives as optional inputs
    if header:
        skiprows=1
    else:
        skiprows=0        

    trans_mat = np.loadtxt(trans_mat_file, delimiter=delimiter, skiprows=skiprows)
    # Extract rotation and transformation matrix from input
    R = trans_mat[0:3,0:3]
    T = trans_mat[0:3,3]
    return R, T

# Move to utils when needed in another function
def insert_columns_with_func(df, labels, subscrpt, func):
    # Insert new columns in dataframe df, next to  columns with labels in list "labels", with name appended with
    # "_subscript", with values given by applying function "func" to the original columns
    if all([label in df.columns for label in labels]):
            # Consider removing the "to_numpy bit and bringing into the function itself"
            values = df[labels].to_numpy()
            # Perform transformation using function
            values_trans = func(values)
            # Insert next to original columns in dataframe
            for i, label in enumerate(labels):
                in_loc = df.columns.get_loc(label)
                df.insert(loc=(in_loc+3), column=(label+subscrpt), value = values_trans[:,i])

if __name__ == "__main__":
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Data_subtracted_disp", out_sub_folder="Data_rotated")  
    # Define file where transformation matrix is stored
    trans_mat_file = os.path.join(Files.parent_path, "image_0000_transformation_matrix.txt")
    # Read in transformation matrix
    R, T = transmat_from_file(trans_mat_file)
    # Otherwise transformation matrices may be defined directly
    # R = np.eye(3)
    # T = np.zeros(3)
    transform_coords(Files, R, T)