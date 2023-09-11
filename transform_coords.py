from operator import is_
import os
import pandas as pd
import numpy as np
import warnings
from FileSeries import *

def transform_coords(Files, label_list, is_displacement = [], R=[], T=[], subscript = "rot"):
    # Read in DIC data from csv files, then apply rotations and translations
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # labels is a list of coordinate labels, or a list of lists of labels if multiple transformation are required
    # is displacement is a list of booleans indicating whether the columns indexed by each entry of label_list
    # are displacements, in which case only a rotation is performed
    # R is a dxd rotation matrix where d is the number of columns indexed by each entry of label_list
    # T is a d-dimensional vector used to translate the coordinates, where d is the number of columns indexed by each entry of label_list
    
    # If a single string has been provided covert to list
    if type(label_list) is not list:
        label_list = list(label_list)
        warnings.warn("Single label has been converted to list as this is the required format")
    # If only list of coordinates has been provided, nest within a further list for iterating
    if type(label_list[0]) is not list:
        label_list = [(label_list)]
    
    if len(is_displacement) == 0:
        is_displacement = [False for labels in label_list]
    elif len(is_displacement) != len(label_list):
        raise Exception("Number of items in is_displacement must match the number of lists in labels")

    # Read data from csv files
    Files.read_data()
    # Apply each required transformation
    for i, labels in enumerate(label_list):
        if is_displacement[i]:
            # If the quantity is a displacement, no need to perform a translation 
            Files.apply_func_to_data(lambda x:rotate_translate(x, coord_label=labels, R=R), labels, subscript, message = "Transforming " + ", ".join(labels))
        else:
            Files.apply_func_to_data(lambda x:rotate_translate(x, coord_label=labels, R=R,T=T), labels, subscript, message = "Transforming " + ", ".join(labels))
    Files.dump()

def rotate_translate(data, coord_label = [], R=[], T=[]):
    # apply a coordinate rotation and transformation to array of coordinates, 
    # data can either be a dataframe with coord_label indicating the columns which represent displacement
    # or the coordinates only, R is a d x d rotation matrix, and T is a  d-dimensional translation vector. 
    # If either R or T are not passed the corresponding transformation is not performed
    # If coord label is not past, coordinates are given by the entire dataframe
    if coord_label:
        # Extract coordinates from the dataframe if a label has been provided, otherwise use the entire
        # dataframe
        xyz = data[coord_label]
    else:
        xyz = data
    xyz = xyz.to_numpy() # Convert to numpy
    if len(R)!=0:
        xyz = (R@xyz.T).T
    if len(T)!=0:
        xyz = xyz + T

    return(xyz)

def transmat_from_file(transmat_file, delimiter=" ", header = False):
    # Load in transformation matrix from a text file
    # trans_mat_file = file where transformation matrix is defined. This is in the format outputted by CloudCompare
    # where the first 3 columns x 3 rows is a rotation matrix, and the 4th column with 3 rows is the translation vector
    # Standard format is delimited by spaces with no header, but can pass alternatives as optional inputs
    if header:
        skiprows=1
    else:
        skiprows=0        

    trans_mat = np.loadtxt(transmat_file, delimiter=delimiter, skiprows=skiprows)
    # Extract rotation and transformation matrix from input
    R = trans_mat[0:3,0:3]
    T = trans_mat[0:3,3]
    return R, T

if __name__ == "__main__":
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Data_subtracted_disp", out_sub_folder="Data_rotated_2")
    # Define file where transformation matrix is stored
    transmat_file = os.path.join(Files.parent_path, "image_0000_transformation_matrix.txt")
    # Read in transformation matrix
    R, T = transmat_from_file(transmat_file)
    labels = [["x", "y", "z"],["x_0","y_0","z_0"],["u","v","w"]]
    is_displacement = [False, False, True]
    transform_coords(Files, labels, is_displacement = is_displacement, R=R, T=T)