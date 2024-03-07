# Gets rid of untrustworthy DIC points outside of the filed of view (around the right corner) based upon 
# Projected coordinates from  get_nat_coords_grid.py

import pandas as pd
import numpy as np
import os
import operator
import re
from FileSeries import *

def trim_data_cond(data, cond, data_del = []):
    # Remove data by row according to logical bitwise cond, and store deleted rows in data_del.
    # If data_del is also provided deleted data is appended to this dataframe
    if len(data_del) == 0:
        # Initialise dataframe for deleted columns if this has not been provided
        data_del = pd.DataFrame(columns=data.columns)

    data_del = pd.concat([data_del, data[cond]], axis=0)
    data = data[~cond]

    return (data, data_del)

def make_na_cond(data, cond, data_del = []):
    # Rather than deleting data from a set, makes NA but retains in the original dataframe
    data[cond] = np.nan
    return(data, data_del)

def str_2_logical(in_string):
    # Convert string containing logical (or arithmetic) expression to evaluate the expression itself
    mapping = {"<"  : operator.lt,
               "<=" : operator.le,
               ">"  : operator.gt,
               ">=" : operator.ge,
               "&"  : operator.and_,
               "|"  : operator.or_,
               "+"  : operator.add,
               "-"  : operator.sub,
               "*"  : operator.mul,
               "/"  : operator.truediv}
    
    return(mapping[in_string])

def partition_by_op(in_string, ops):
    # Partition string according to specified operators
    out_list = [in_string]
    out_ops = []
    for op in ops:
        for i, string in enumerate(out_list):
            if op in string:
                # I've already added spaces around the minus sign to distinguish between minus and subtract
                if not op == " - ": 
                    # Pad with spaces if not already
                    out_list[i] = string.split(op.center(len(op)+2))
                else:
                    out_list[i] = string.split(op)
                out_list = [item for item in out_list]
                # Flatten the list to avoid nesting
                flt_list = []
                for item in out_list:
                    if isinstance(item, list):
                        flt_list.extend(item)
                    else:
                        flt_list.append(item)
                out_list = flt_list

    if len(out_list) == 1:
        out_list = out_list[0]
    
    # Use regular expressions to return a list of the operators in the order in which they appear
    for char in ["|","+","-","*"]:
        ops = [op.replace(char,"\\" + char) for op in ops]

    pattern = r"|".join(ops)
    out_ops = re.findall(pattern, in_string)
    return(out_list, out_ops)

def get_log_index(data, in_string):
    # Get logical index for applying a set of instructions, contained in a string, to a dataframe
    # Avoid using "eval" as this is bad practice for security reasons
    # Partition by and and or operations
    indiv_conds, linking_ops = partition_by_op(in_string, ["&","|"])
    # how deep is indiv_ops? If only one convert to list for looping
    if len(indiv_conds[0]) == 1:
        indiv_conds = [indiv_conds]

    indiv_inds = []
    for cond in indiv_conds:
        # Now partition string by inequalities
        cond, indiv_op = partition_by_op(cond, ["<=",">=",">","<"]) # Check less than or equal to first otherwise they get wrongly identified
        # If a simple condition the second term should be a float. If this doesn't work further partitioning is required
        try:
            cond[1] = float(cond[1])
            # Append logical index to list
            indiv_inds.append(str_2_logical(indiv_op[0])(data[cond[0]],cond[1]))
        except:
            # Further partition second term for arithmetic operations
            # Add spaces to " - " to distinguish subtract from minus
            complex_cond, complex_cond_ops = partition_by_op(cond[1], ["+"," - ","*","/"])
            # Extract the data from thte dataframe and apply arithmetic operations,
            # then compute the logical index
            data_cond = data[complex_cond[0]]
            numeric_vals = [float(item) for item in complex_cond[1:]]
            for op, val in zip(complex_cond_ops, numeric_vals):
                data_cond = str_2_logical(op)(data_cond,val)
            indiv_inds.append(str_2_logical(indiv_op[0])(data[cond[0]],data_cond))
    
    # Get function for linking and/or term
    linking_ops = [str_2_logical(op) for op in linking_ops]
    # Link together logical indices with and/or as required
    out_ind = indiv_inds[0]
    if len(linking_ops) > 0:
        for op, link in zip(indiv_inds[1:],linking_ops):
            out_ind = link(out_ind, op)

    return(out_ind)

def trim_bounds(Files, bound_list, dropna = False, out_path_del = []):    
    # Main function for trimmming data
    # Use a different function depending on if nas are to be retained
    if dropna:
        filter_func = trim_data_cond
    else:
        filter_func = make_na_cond

    for file in Files.files:
        for bound in bound_list:
            file.data, data_del = filter_func(file.data, get_log_index(file.data, bound))

        # Write output to csvs
        print("writing " + file.dst)
        file.data.to_csv(file.dst, sep=",", index=True)
        # Need to implement capability to store deleted data
        #if dropna and out_path_del:
            #if out_path_del not in os.listdir():
            #    os.mkdir(out_path_del)
            #data_del.to_csv(os.path.join(wd, out_folder_del, file), sep=",", index=False)
            #data_del.to_csv(out_path_del, sep=",", index=False)

if __name__ == "__main__":
    folder = "..\\CS02P\\DIC\\Left_Camera_Pair"
    # IMPORTANT ISSUE: MIGHT NOT BE A PROBLEM WITH THIS CODE, PROBABLY MORE AN ISSUE WITH FILESERIES
    # BUT WHEN I RUN THIS AND INTEGERS ARE READ IN FROM CSVS THEY ARE INTERPRETED AS FLOATS
    # SORT IT OUT
    Files = FileSeries(folder=folder,in_sub_folder="Working_Folder", out_sub_folder="Trimmed_Rad")
    Files.read_data()
    out_path_del = os.path.join(folder, "Rad_trimmed_del")
    # coord_labels = ["x_proj","y_proj","z_proj"]
    bound_list = ["y_proj < 5.0", # points in the obscured radius of the ends
                  "z_proj >= 185.0 & z_proj <= 235.0 & y_proj < 11.25", # Points in the obscured radius in the central region
                  "z_proj >= 60.0 & z_proj < 185.0 & y_proj < z_proj / 20.0 + 2.0",
                  "z_proj > 235.0 & z_proj <= 360.0 & y_proj < z_proj / -20.0 + 23.0"]
                  
    trim_bounds(Files, bound_list, out_path_del=out_path_del)