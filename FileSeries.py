import os
import pandas as pd
import warnings
from natsort import natsorted
import numpy as np
from utils import *
from str2bool import *

class FileSeries:
    # Contains properties of a series of .csv files associated with a particular test
    def __init__(self,folder="", in_sub_folder="", out_sub_folder="", del_sub_folder=""):
        # Perhaps pass folder as a kwargs? If not passed assume that data is in the working directory
        self.folder = folder # Parent folder in which the data is stored
        self.in_sub_folder = in_sub_folder # Subfolder in which the data is stored
        if out_sub_folder == "":
            self.out_sub_folder = self.in_sub_folder
        else:
            self.out_sub_folder = out_sub_folder
        self.del_sub_folder = del_sub_folder # Subfolder where deleted data will be written, if reqired
        self.get_paths()
        # Create a directory for writing output if required
        make_out_folder(self.out_path)
        # Create a list of all files in the in_sub_folder
        self.get_files()

    def down_sam(self, rate, rename_files = False):
        # Only read data from a subset of files
        self.files = [file for i, file in enumerate(self.files) if i%rate == 0]
        # Need to rename the files after?
        if rename_files:
            for i, file in enumerate(self.files):
                file.out_filename = file.out_filename.rstrip("0123456789"+file.ext)+str(i)+".csv"
                file.update_dst()
    
    def trunctate_data(self, start = 0, end = -1):
        # Trim the data to a specific range
        self.files = [self.files[i] for i in range(start,end)]

    def get_paths(self):
        # Get full path to which data is contained
        wd = os.path.abspath("") # Get working directory
        self.parent_path = os.path.join(wd,self.folder)
        self.in_path = os.path.join(wd,self.folder,self.in_sub_folder)
        self.out_path = os.path.join(wd,self.folder,self.out_sub_folder)
        if self.del_sub_folder != "":
            self.del_path = os.path.join(wd,self.folder,self.del_sub_folder)

    def get_files(self):
        # Create a list of File objects
        if hasattr(self,"del_path"):
            self.files = [File(filename, self.in_path, self.out_path, del_path=self.del_path) for filename in natsorted(os.listdir((self.in_path)))]
        else:
            self.files = [File(filename, self.in_path, self.out_path) for filename in natsorted(os.listdir((self.in_path)))]

    def read_data(self, dropna = False, **kwargs):
        # Consider adding a progress bar
        # use kwargs to pass keywords to pandas read_csv function
        print("Reading data")
        for i, File in enumerate(self.files):
            print(i)
            File.read_file(dropna, **kwargs)

    def extract_qoi(self, qoi, new_names = [], dropna = True, **kwargs):
        # Extracts columns of interest from a csv files, then removes rows with missing values
        # QoI = list of strings indicating which columns are to be extracted
        # new_names = list of updated names which are to be applied to the columns. Must be the same length as QoI
        print("Filtering data")
        for i, File in enumerate(self.files):
            print(i)
            File.filter_data(qoi, new_names, dropna, **kwargs)

    def filter_by_cond(self, conds, drop_cond = False, dropna = False, **kwargs):
        # Filters the data using some logical index given in a string
        # cond = list of strings containing logical conditions
        # drop_cond = if true, remove the data which satisfies the condition, otherwise default to keep
        for cond in conds:
            print("Filtering data by: " + cond)
            for i, File in enumerate(self.files):
                # get logical index for the condition
                print(i)
                bool_ind = get_log_index(File.data, cond)
                if drop_cond:
                    bool_ind = ~bool_ind
            
                if dropna:
                    File.trim_data_cond(bool_ind, **kwargs)    
                else:
                    File.make_na_cond(bool_ind, **kwargs)

    def update_datatype(self, dtype, columns = [], index = []):
        # Update datatpye of columns specified in "columns" or "index" lissts, to type specified by 
        # string "dtype".  Columns may be a list of column names, index list of interger or indices.
        # If empty is applied to whole dataframe
        print("Converting datatypes")
        for i, File in enumerate(self.files):
            print(i)
            File.update_datatype(dtype, columns = columns, index = index)

    def apply_func_to_data(self, func, labels = [], in_sub = [], out_sub = "", rel_pos = 1, out_cols = [], message = [],  insert_by_label = True, first_file_only = False):
        # Applies a function to every dataframe, and inserts new columns with the results. 
        # Primarily used for coordinate transformations
        # func = function which is applied to the data given as lambda function
        # labels = labels of columns to which function is applied
        # subscript = subscript added to the label of the new columns
        # rel_pos = relative position (to the right of the input column) at which new column is inserted
        # out_cols = list of integer indices of columns to which data will be added if not by label. Default value is to end of dataframe
        # message = optional message printed to the terminal
        # insert_by_label = Boolean indicating whether to insert output columns relative to input columns (True) or by index (False)
        # first_file_only = Boolean indicating to apply function only to the first file, and duplicate output for all other files. Use if output expected to be identical.
        if message:
            print(message)
        if first_file_only:
            # Apply function to the first file, and duplicate output to all other files
            print("0")
            values = func(self.files[0].data)
            if insert_by_label:
                if not labels:
                    raise Exception("To specify new column position relative to label, please input list of labels as labels=")
                for i, File in enumerate(self.files):
                    if i > 0:
                        print(i)
                    if values.shape[0] != File.data.shape[0]:
                        raise Exception("Datasets are different sizes across files, cannot duplicate output from first file across all datasets")
                    File.insert_col_by_label(values, labels, in_sub = in_sub, out_sub = out_sub, rel_pos = rel_pos)
            else:  
                if not out_cols: 
                    raise Exception("To specify new column position with index, please return list of indices as out_cols=")             
                for i, File in enumerate(self.files):
                    if i > 0:
                        print(i)
                    if values.shape[0] != File.data.shape[0]:
                        raise Exception("Datasets are different sizes across files, cannot duplicate output from first file across all datasets")
                    File.insert_col_by_loc(values, out_cols=out_cols)

        else:
            for i, File in enumerate(self.files):
                print(i)
                if insert_by_label:
                    if not labels:
                        raise Exception("To specify new column position relative to label, please input list of labels as labels=")
                    File.insert_col_with_func_by_label(func, labels, in_sub = in_sub, out_sub = out_sub, rel_pos = rel_pos)
                else:
                    if not out_cols:
                        raise Exception("To specify new column position with index, please return list of indices as out_cols=")
                    File.insert_col_with_func_by_loc(func, out_cols=out_cols)
                # Possibly also add option to just replace the column if specified by Boolean

    def dump(self, dropna = False, **kwargs):
        # For just renaming it's quicker to copy. Incorporate this here when ready...
        print("Writing Processed data")
        for i, File in enumerate(self.files):
            print(i)
            File.write_data(dropna = dropna, **kwargs)

    def dump_data_del(self, **kwargs):
        # Use kwargs to pass arguments to pandas to_csv()
        print("Writing deleted data")
        if not hasattr(self, "del_path"):
            self.del_path = self.out_path+"_del"
        make_out_folder(self.del_path) # Create subfolder if it doesn't already exist

        for i, File in enumerate(self.files):
            print(i)
            # Update file definition if a destination for deleted data hasn't been given
            if not hasattr(File, "del_dst"):
                File.update_del_dst(self.del_path) 
            File.write_del_data(**kwargs)

    def dump_nas(self, **kwargs):
        # Use kwargs to pass arguments to pandas to_csv()
        # Copy rows with nas to data_del
        for File in self.files:
            File.flag_nas_as_del()
        self.dump_data_del(**kwargs)

class File:
    # Properties of an individual file within a FileSeries
    def __init__(self, filename, in_path, out_path, del_path = ""):
        self.in_filename = filename
        self.out_filename = filename # Default output filename is the same as the input
        self.in_path = in_path
        self.ext = os.path.splitext(filename)[-1] # Get file extension
        self.out_path = out_path
        self.src = os.path.join(self.in_path, filename)
        if self.ext == ".vtk":
            self.out_filename = os.path.splitext(filename)[0] + ".csv"
            #self.dst = os.path.join(self.out_path, self.out_filename)
        #else:
        #    self.dst = os.path.join(self.out_path, filename)
        self.dst = os.path.join(self.out_path, self.out_filename)
        if del_path != "":
            self.del_path = del_path
            self.del_dst = os.path.join(self.del_path, self.out_filename)    
        self.data = [] # add data from csv files

    def update_dst(self):
        # Update the destination address if this is renamed
        self.dst = os.path.join(self.out_path, self.out_filename)

    def update_del_dst(self, del_path):
        # Update the destination address if this is renamed
        self.del_path = del_path
        self.del_dst = os.path.join(self.del_path, self.out_filename)

    def read_file(self, dropna, **kwargs):
        # Use kwargs to pass options to pandas read_csv function
        if self.ext == ".vtk":
            self.data = vtk_to_pandas(self.src)
        else:
            self.data = pd.read_csv(self.src, **kwargs)
        if dropna:
            self.data.dropna()
        self.n_points = self.data.shape[0]

    def filter_data(self, qoi, new_names = [], dropna = True, **kwargs):
        self.data = self.data.filter(items=qoi, **kwargs)
        # Rename columns if required
        if new_names:
            if len(self.data.columns) == len(new_names):
                self.data.columns = new_names 
            else:
                warnings.warn("Specified number of new_names did not match the selected number of Colummns. Could not rename.")
        # Drop na values if required
        if dropna:
            self.data = self.data.dropna(ignore_index=True)
            
        self.n_points = self.data.shape[0]

    def make_na_cond(self, cond, na_col = []):
        # Make data points na based upon some logical index, which must have the same number of rows as self.data. 
        # Used to mark datapoints for later deletion
        # cond is a logical index
        # na_col is column used to flag nas. If not specified all columns are set to na
        if na_col:
            self.data.loc[~cond, na_col] = np.nan
        else:
            self.data[~cond] = np.nan

    def trim_data_cond(self, cond, store_data_del = False):
        # Remove data by row based upon some logical index, which must have the same number of rows as self.data
        # If specified, retain deleted data 
        if store_data_del:
            if not hasattr(self, "data_del"):
                self.data_del = pd.DataFrame(columns=self.data.columns)
            self.data_del = pd.concat([self.data_del, self.data[~cond]], axis=0)
        self.data = self.data[cond]

    def flag_nas_as_del(self):
        # Create a data_del attribute containing all rows of self.data containing na values, to track
        # points which are deleted
        if not hasattr(self, "data_del"):
            self.data_del = pd.DataFrame(columns=self.data.columns)        
        self.data_del = pd.concat([self.data_del, self.data[self.data.isna().any(axis=1)]])

    def update_datatype(self, dtype, columns = [], index = []):
        # Update datatpye of columns specified in "columns" or "index" lissts, to type specified by 
        # string "dtype".  Columns may be a list of column names, index list of interger or indices.
        # If empty is applied to whole dataframe        
        if columns:
            self.data[columns].astype(dtype)
        elif index:
            self.data.iloc[:,index].astype(dtype)
        else:
            self.data.astype(dtype)

    def insert_col_by_label(self, values, labels, in_sub = [], out_sub = "", rel_pos = 1):
        # Insert output to the right of input columns with labels in list "labels" + _"in_sub", with 
        # name appended by "_out_sub", using values spefieid in "values".
        # optional input rel_pos gives the relative position of the new column compared to the original
        # Only works if function output has same number of columns as input
        if in_sub:
            in_labels = ["_".join((label, in_sub)) for label in labels]
        else:
            in_labels = labels
        out_labels = ["_".join((label, out_sub)) for label in labels]
        if all([label in self.data.columns for label in in_labels]):
            if values.shape[1] != len(in_labels):
                raise Exception("The number of columns outputted by the function do not match the specified number of labels")
            
            # Insert transformed data into the dataframe
            for i, (in_label, out_label) in enumerate(zip(in_labels, out_labels)):
                in_loc = self.data.columns.get_loc(in_label)
                self.data.insert(loc=(in_loc+rel_pos), column=out_label, value = values[:,i])
        else:
            raise Exception("The specified columns could not be found in the data")

    def insert_col_with_func_by_label(self, func, labels, **kwargs):#in_sub = [], out_sub = "", rel_pos = 1):
        # given by applying function "func" to the original columns
        # Insert output to the right of input columns with labels in list "labels" + _"in_sub", with 
        # name appended by "_out_sub", and values given by applying function "func" to the original columns
        # optional input rel_pos gives the relative position of the new column compared to the original
        # Only works if function output has same number of columns as input
        
        # Apply the function to the data
        values = func(self.data)        
        self.insert_col_by_label(values, labels, **kwargs)            
        
    def insert_col_by_loc(self, values, out_cols = []):
        if len(out_cols) < values.shape[1]:
            # If insufficient output column indices are specified then place at end of Dataframe
            out_cols.extend([-1 for i in range(values.shape[1]-len(out_cols))])
        elif len(out_cols) > values.shape[1]:
            warnings.warn("More output column locations specified than are outputted by the function, the last {} values are ignored".format(len(out_cols)-values.shape[1]))
        for loc, column in zip(out_cols, values):
            if loc == -1:
                # If needed at end can just create new column rather than using insert
                self.data[column] = values[column]
            else:
                self.data.insert(loc=loc, column=column, value = values[column])
    
    def insert_col_with_func_by_loc(self, func, out_cols = []):
        # Run function
        values = func(self.data)
        self.insert_col_by_loc(values, out_cols)

    def write_data(self, sep = ",", index = False, dropna = False, **kwargs):
        # Use kwargs to pass writing options to pandas
        if dropna:
            self.data.dropna().to_csv(self.dst, sep = sep, index = index, **kwargs)
        else:
            self.data.to_csv(self.dst, sep = sep, index = index, **kwargs)

    def write_del_data(self, sep = ",", index = False, **kwargs):
        # Write deleted data to folder
        if not hasattr(self,"data_del"):
            raise Exception("No deleted data exists. If tracking with NAs, use dump_nas() instead")
        self.data_del.to_csv(self.del_dst, sep = sep, index = index, **kwargs)