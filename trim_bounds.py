# Trims the data to remove a set of points from outside the region of interest
# Trims based upon a set of logical condions past via a list of strings. Within
# each string indices, operators and numerical values must be separated by a 
# single space

from FileSeries import *

def trim_bounds(Files, bound_list, drop_cond = False, dropna = False, **kwargs):
    # Main function for trimmming data
    Files.filter_by_cond(bound_list, drop_cond = drop_cond, dropna = dropna, **kwargs)

if __name__ == "__main__":
    # Select subset of data within specified bounds

    # Create FileSeries object and read in data
    folder = "..\\CS02P\\DIC\\Right_Camera_Pair"
    Files = FileSeries(folder=folder,in_sub_folder="Processed_Data_Working", out_sub_folder="Trimmed_Rad", del_sub_folder = "Del_Data")
    Files.read_data()

    # Define list of bounds defining where to trim data
    #bound_list = ["y_proj < 5.0", # points in the far radius of the ends
    #              "z_proj >= 185.0 & z_proj <= 235.0 & y_proj < 11.25", # Points in the obscured radius in the central region
    #              "z_proj >= 60.0 & z_proj < 185.0 & y_proj < z_proj / 20.0 + 2.0",
    #              "z_proj > 235.0 & z_proj <= 360.0 & y_proj < z_proj / -20.0 + 23.0"]
    #bound_list = ["y_proj > 145", #points in the far radius of the ends
    #              "z_proj >= 185.0 & z_proj <= 235.0 & y_proj > 138.75", # Points in the far radius in the central region
    #              "z_proj >= 60.0 & z_proj < 185.0 & y_proj > -z_proj / 20.0 + 148.0", # Point in the far radius in the ramp
    #              "z_proj > 235.0 & z_proj <= 360.0 & y_proj > z_proj / 20.0 + 127.0"]

    # bound_list = ["z_proj >= 407.5 & z_proj >= 407.5], #points in the far radius of the ends
    #              "z_proj >= 185.0 & z_proj <= 235.0 & y_proj > 138.75", # Points in the far radius in the central region
    #              "z_proj >= 60.0 & z_proj < 185.0 & y_proj > -z_proj / 20.0 + 148.0", # Point in the far radius in the ramp
    #              "z_proj > 235.0 & z_proj <= 360.0 & y_proj > z_proj / 20.0 + 127.0"]

    trim_bounds(Files, bound_list, drop_cond=True, dropna=False, na_col = "Time")
    # Write output
    Files.dump_nas() # Write deleted data to separate folder, tracked using nas
    # Files.dump_data_del() Alternative for if deleting data on the fly rather than tracking with nas
    Files.dump(dropna=True) 