import os
import pandas as pd
import numpy as np
from FileSeries import *
from utils import model_force_from_json
from downsample_data import * # Won't need this in long run, put in manage workflow

# Interpolate data to specified applied load values
# In future consider adding other methods, e.g. moving average?
def interp_data(Files, model_force, interp_by_col = "Load", downsam = False, downsam_rate = 1, dropna = True, output_all = True, output_frames = True):
    # Extract forces from each image
    force_data = [-file.data[interp_by_col][0] for file in Files.files]

    # Loop over each requested applied force and find the image with force either side of this value
    before_ind = []
    after_ind = []
    for force_inc in model_force:
        # Check if there is experimental data lower than the reqested applied force
        if any([force < force_inc for force in force_data]):
            # Append the last entry for which the force is lower than the requested value
            before_ind.append([i for i, force in enumerate(force_data) if force < force_inc][-1])
        else:
            # Otherwise just use the first datapoint
            # before_ind.append(0)
            # before_ind.append([])
            # before_ind.append(np.NaN)
            before_ind.append(-1) # replace later with nas, but after converting to integer
            # Consider allowing for empty value
        
        # Append the first entry for which the experimental force is higher than the requested value to the list of "after" indices
        after_ind.append([i for i, force in enumerate(force_data) if force >= force_inc][0])
    
    interp_frame = pd.DataFrame({"Model_Force" : model_force})
    interp_frame["Before_Force"] = pd.Series([force_data[ind] if ind != -1 else pd.NA for ind in before_ind])
    interp_frame["After_Force"] = pd.Series([force_data[ind] for ind in after_ind])
    interp_frame["Before_Image"] = pd.Series(before_ind,dtype=int)
    interp_frame["Before_Image"].replace(-1, pd.NA, inplace=True)
    interp_frame["After_Image"] = pd.Series(after_ind,dtype=int)
    interp_frame["Ratio"] = (interp_frame["Model_Force"] - interp_frame["Before_Force"])/(interp_frame["After_Force"] - interp_frame["Before_Force"])
    # If there is no before then this will give an infinite value, as denominator is zero
    # Just replace with 1 to use the after image
    interp_frame.replace([np.inf, -np.inf], 1.0, inplace = True)
    # Identify the data corresponding to the image before and after each force, and linearly
    # interpolate. For this to work the index of the data should be retained throughout the processing
    # otherwise it will be impossible to track points from one image to the next across load steps
    for i, row in enumerate(interp_frame.itertuples()):
        # Possibility of using index from dataframe rather than enumerate (provided this works, which it may not)
        # before_DIC = pd.read_csv(os.path.join(image_folder,row["Before_Image"]))
        # if not np.isnan(row.Before_Image):
        if not pd.isna(row.Before_Image):
            before_DIC = Files.files[row.Before_Image].data
            after_DIC = Files.files[row.After_Image].data
            # after_DIC = pd.read_csv(os.path.join(image_folder,row["After_Image"]))
            if before_DIC.shape[0] != after_DIC.shape[0]:
                raise Exception("Different number of datapoints across images, cannot interpolate")

            # Interpolate between the two data points
            intp_DIC = (after_DIC - before_DIC)*row.Ratio + before_DIC
            # Don't interpolate integers
            # Ideally I'd detect integers automatically, but a previous bit of code has 
            # written these to csv as float. Hopefully I'll fix this in another iteration,
            # at which point I can update the code. For now do it manually.
            # print(before_DIC.dtypes)
            # print(before_DIC.dtypes != "float64")
            #intp_DIC[before_DIC.dtpyes != "float64"] = before_DIC[before_DIC.dtpyes != "float64"]
            intp_DIC[["Element","Conv_Iteration"]] = before_DIC[["Element","Conv_Iteration"]]
            intp_DIC["Increment"] = i + 1
            intp_DIC["Compressive Force"] = row.Model_Force
            if downsam:
                print("this shouldn't run")
                intp_DIC = downsample_data(intp_DIC, sort = True, sort_col = sort_col, method = "uniform", rate = downsam_rate)
            # Later this will be incorporated into the FileSeries object, but would need separate entry for
            # interpolated data. Keep here for now
            if dropna:
                intp_DIC = intp_DIC.dropna()
            # Output individual frames if required
            if output_frames:
                intp_DIC.to_csv(os.path.join(Files.files[0].out_path, "Image_Inc_" + str(i) + ".csv"),sep=",",index=True)
            # If required append to dataframe containing all interpolated output
            if output_all:
                try:
                    all_frames_DIC = pd.concat((all_frames_DIC,intp_DIC.dropna()),axis=0)
                except:
                    all_frames_DIC = pd.DataFrame(columns=intp_DIC.columns.values)
                    all_frames_DIC = pd.concat((all_frames_DIC,intp_DIC.dropna()),axis=0)

    if output_all:
        # Write to csv
        all_frames_DIC.to_csv(os.path.join(Files.in_path, "Interpolated_DIC.csv"), sep=",", index=True)

if __name__ == "__main__":
    folder = "..\\CS02P\\DIC\\Left_Camera_Pair"
    Files = FileSeries(folder=folder,in_sub_folder="Trimmed_Rad", out_sub_folder="Interpolated_Data")
    Files.read_data()
    # truncate the data to before failure occurs
    Files.trunctate_data(end = 336)
    interp_by_col = "Load"
    model_force = [15*i for i in range(15)]
    # alternative code where force is extracted from the model output data
    model_out_json = "C:\\Users\\cs2361\\Documents\\CSpar_Calibration\\inputs\\LHSDesign50x5_output_struct_200kN.json"
    model_force = model_force_from_json(model_out_json)
    print(model_force)
    downsam = False
    downsam_rate = 8
    print(model_force)
    interp_data(Files, model_force, interp_by_col = "Load", downsam = downsam)
    sort_col = ["Element","h"]