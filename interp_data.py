import os
import pandas as pd
import numpy as np
import json

# Interpolate data to specified applied load values
# In future consider adding other methods, e.g. moving average?
def downsample_data(data, sort = False, sort_col = ["Element"], method = "uniform", rate = 2):
    # Downsample the data to reduce sample size
    n_data = data.shape[0]
    # Sort into order of ascending element number such that points in close proximity are grouped
    if sort:
        data.sort_values(sort_col,inplace=True)
    if method == "uniform":
        data = data.iloc[range(0,n_data,rate)]
    elif method == "random":
        data = data.sample(frac = 1.0/rate)
    return(data)

wd = os.getcwd()
folder = "..\\..\\..\\..\\Load Displacement Curve"
file = "Image.csv"
image_folder = "Data_rad_trimmed\\"
down_sam = True # Downsample to thin data?
downsam_rate = 8

data = pd.read_csv(os.path.join(wd,folder,file), sep=";")
data = data.filter(items = ["File", " Force [kV]", " Displacement [mm]"])
data.columns = ["Image", "Force", "Displacement"]
images = ["Image_" + file.partition("Image_")[2].replace("_0.tiff",".csv") for file in data["Image"]]
data["Image"] = pd.Series(images)
data["Compressive_Force"] = -data["Force"]
# Not sure, how to play this what with Geir subtracting the first displacement component. Consider cheating by
# subtracting the first force to get zero - this is cheating but maybe ok...
data["Zeroed_Force"] = data["Compressive_Force"] - data["Compressive_Force"][0]
# Could try both

data.to_csv("Force-Displacement.csv", sep=",",index=False)

# Manually truncate the data to before failure occurs
data = data.filter(items = range(181), axis = 0)

# Interpolate DIC to fixed load increments to match model output
# model_output_location = "C:\\Users\\cs2361\\Documents\\CSpar_Calibration\\inputs\\LHSDesign40x4_output_struct.json"
model_output_location = "C:\\Users\\cs2361\\Documents\\CSpar_Calibration\\inputs\\LHSDesign70x7_downsam_2.json"
with open(model_output_location, "r") as f:
    # Load in string from file
    model_dict = json.loads(f.readline())

model_load = [frame["RFs"][2] for frame in model_dict["Sample"][0]["Frame"]]

# Consider keeping filenames the same for automation, rather than stripping the word "Image" out
before_ind = []
after_ind = []
# Only extract DIC for the third last frame
# model_load = [model_load[-3]]
#print(model_load)
for load_inc in model_load:
    # We want to extract the image either side of the current load increment then interpolate between the 2, if possible
    if not data[data["Compressive_Force"]<load_inc].empty: 
#        print(data[data["Compressive_Force"]<load_inc].index[-1])
        before_ind.append(data[data["Compressive_Force"]<load_inc].index[-1])
        # print(data[data["Compressive_Force"]<load_inc].iloc[-1]["Compressive_Force"])
    else:
        before_ind.append(0)
        #print("First entry")

    after_ind.append(data[data["Compressive_Force"]>=load_inc].index[0])
    # print(data[data["Compressive_Force"]>=load_inc].iloc[0]["Compressive_Force"])
        
interp_data = pd.DataFrame({"Model_Force" : model_load})
interp_data["Before_Force"] = data["Compressive_Force"].filter(items = before_ind, axis = 0).values
interp_data["After_Force"] = data["Compressive_Force"].filter(items = after_ind, axis = 0).values
interp_data["Before_Image"] = data["Image"].filter(items = before_ind, axis = 0).values
interp_data["After_Image"] = data["Image"].filter(items = after_ind, axis = 0).values
interp_data["Ratio"] = (interp_data["Model_Force"] - interp_data["Before_Force"])/(interp_data["After_Force"] - interp_data["Before_Force"])
# If there is no before then this will give an infinite value, as denominator is zero
# Just replace with 1 to use the after image
interp_data.replace([np.inf, -np.inf], 1.0, inplace = True)

print(interp_data)

# out_folder = "Interpolated_Data"
out_folder = "adasd"
if out_folder not in os.listdir(os.getcwd()):
    os.mkdir(out_folder)

# For this to work, I need to have retained the index of the data throughout the processing
# Otherwise I won't be able to track points across load steps
for i,row in interp_data.iterrows():
    print(i)
    before_DIC = pd.read_csv(os.path.join(image_folder,row["Before_Image"]))
    after_DIC = pd.read_csv(os.path.join(image_folder,row["After_Image"]))
    if before_DIC.shape[0] != after_DIC.shape[0]:
        print(i)
        print(before_DIC)
        print(after_DIC)
        asdsad

    # Interpolate between the two data points
    intp_DIC = (after_DIC - before_DIC)*row["Ratio"] + before_DIC
    # Don't interpolate integers...
    intp_DIC[["Element","Conv_Iteration"]] = before_DIC[["Element","Conv_Iteration"]]
    intp_DIC["Increment"] = i + 1
    intp_DIC["Compressive Force"] = row["Model_Force"]
    # Downsample the data if necessary
    if down_sam:
        sort_col = ["Element","h"]
        intp_DIC = downsample_data(intp_DIC, sort = True, sort_col = sort_col, method = "uniform", rate = downsam_rate)

    intp_DIC.dropna().to_csv(os.path.join(out_folder, "Image_Inc_" + str(i) + ".csv"),sep=",",index=True)
    if i == 0:
        all_frames_DIC = pd.DataFrame(columns=intp_DIC.columns.values)

    all_frames_DIC = pd.concat((all_frames_DIC,intp_DIC.dropna()),axis=0)

print(all_frames_DIC)
# Write to csv
all_frames_DIC.to_csv("Interpolated_DIC.csv", sep=",", index=True)