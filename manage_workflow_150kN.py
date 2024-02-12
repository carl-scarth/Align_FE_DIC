# import os
from pydoc import describe
import matplotlib.pyplot as plt
import numpy as np

from FileSeries import *
from rename_files import *
from subtract_displacement import *
from transform_coords import *

def process_loadstep(folder, in_sub_folder, QoI, new_names, skip_images = 0, mean_stats = [], sd_stats = []):
    print(in_sub_folder)
    out_sub_folder = in_sub_folder + "_processed"
    test_data = FileSeries(folder=folder, in_sub_folder=in_sub_folder, out_sub_folder=out_sub_folder)
    [remove_suffix(File) for File in test_data.files]
    test_data.read_data(sep=",")
    test_data.extract_qoi(QoI, new_names=new_names, dropna = False)
    # Clear all rows containing nas as all values are corrupted for these entries
    for file in test_data.files:
        file.data[file.data.isna().any(axis=1)] = np.nan
    # Calculate point-wise statistics across all images
    # Create dictionary for storing output
    point_dict = {column : pd.DataFrame(test_data.files[skip_images].data[column].values.reshape(1,-1)) for column in test_data.files[skip_images].data.columns.values}
    for key in point_dict.keys():
        for file in test_data.files[skip_images+1:]:
            point_dict[key] = pd.concat((point_dict[key], pd.DataFrame(file.data[key].values.reshape(1,-1))), axis=0)

    # Mean and standard deviation across each point (columns rearranged as rows and stacked vertically)
    point_means = pd.DataFrame({key : item.mean() for key, item in point_dict.items()})
    point_sds = pd.DataFrame({key : item.std() for key, item in point_dict.items()})
    print("Summary statistics of pointwise mean across images")
    print(point_means.describe())
    print("Summary statistics of pointwise standard deviation across images")
    print(point_sds.describe())
    # Store descriptive statistics
    mean_stats.append(point_means[["u", "v", "w", "uvw_magnitude","EI","EII"]].describe())
    sd_stats.append(point_sds[["u", "v", "w", "uvw_magnitude","EI","EII"]].describe())
    # Write point-wise means and standard deviations to csv
    point_means.to_csv("\\".join([folder,out_sub_folder,"point_means.csv"]), sep=",", index=False)
    # Add mean coordinates for plotting
    point_sds[["x_mean","y_mean","z_mean"]] = point_means[["x","y","z"]]
    point_sds.to_csv("\\".join([folder,out_sub_folder,"point_sds.csv"]), sep=",", index=False)

    # Drop nas and write processed data to csv
    for file in test_data.files:
        file.data.dropna(inplace=True)
    test_data.dump()
    return mean_stats, sd_stats, test_data

QoI = ["coor.X [mm]",
       "coor.Y [mm]",
       "coor.Z [mm]",
       "disp.Horizontal Displacement U [mm]",
       "disp.Vertical Displacement V [mm]",
       "disp.Out-Of-Plane: W [mm]",
       "disp.Displacement Magnitude [mm]",
       "strain.Maximum Principal Strain: EI [ ]",
       "strain.Minimum Principal Strain: EII [ ]",
       "stats.Sigma [mm]"]
    
new_names = ["x", "y", "z", "u", "v", "w", "uvw_magnitude", "EI", "EII","sigma"]

sub_folders = ["0kN", "50kN", "100kN", "150kN"]
folder = "..\\150kN Data\\Camera_Pair_0_3"

# Deal with the first subfolder separately
in_sub_folder = sub_folders.pop(0)
mean_stats, sd_stats, test_data = process_loadstep(folder, in_sub_folder, QoI, new_names, skip_images = 1)

# Plot statistics across each image (only meaningful for 0kN)
means = pd.DataFrame([file.data.dropna().mean() for file in test_data.files[1:]])
sds = pd.DataFrame([file.data.dropna().std() for file in test_data.files[1:]])
skews =  pd.DataFrame([file.data.dropna().skew() for file in test_data.files[1:]])
means[["u","v","w"]].plot(title = "Mean across each image at 0kN", xlabel="Image",ylabel="Mean (mm)")
sds[["u","v","w","uvw_magnitude"]].plot(title = "Standard deviation across each image at 0kN", xlabel="Image",ylabel="Standard deviation (mm)")
skews[["u","v","w"]].plot(title = "Skewness across each image at 0kN", xlabel="Image",ylabel="Skew")
# medians = pd.DataFrame([file.data.median() for file in test_data.files[1:]])
# medians[["u","v","w"]].plot(title = "Median across each image at 0kN", xlabel="Image",ylabel="Median (mm)")
# iqr = pd.DataFrame([file.data.quantile(q=0.75) - file.data.quantile(q=0.25) for file in test_data.files[1:]])
# iqr[["u","v","w","uvw_magnitude"]].plot(title = "Inter-quartile range across each image at 0kN", xlabel="Image",ylabel="Interquartile range (mm)")

print("Summary statistics of mean and standard deviation across each image")
print("Mean mean")
print(means.mean())
print("Standard deviation of mean")
print(means.std())
print("Mean of tandard deviation")
print(sds.mean())
print("Maximum standard deviation")
print(sds.max())

# Concatenate data into a single frame then calculate mean and standard deviation
for i, file in enumerate(test_data.files[1:]):
    if i == 0:
        all_data = file.data
    else:
        all_data = pd.concat((all_data,file.data),axis=0)

print("Mean across all data points")
print(all_data.mean())
print("Standard deviation across all data points")
print(all_data.std())

# Repeat relevant steps for all loads and produce a summary plot
for in_sub_folder in sub_folders:
    mean_stats, sd_stats, *_ = process_loadstep(folder, in_sub_folder, QoI, new_names=new_names, mean_stats = mean_stats, sd_stats=sd_stats)

# Plot key statistics at different loads
x_plot = [0, 50, 100, 150]
for i in range(len(x_plot)):
    if i == 0:
        mean_frame = pd.DataFrame(mean_stats[i].loc["mean"].values.reshape(1,-1), columns=mean_stats[i].columns.values)
        meanstd_frame = pd.DataFrame(mean_stats[i].loc["std"].values.reshape(1,-1), columns=mean_stats[i].columns.values)
        sd_frame = pd.DataFrame(sd_stats[i].loc["mean"].values.reshape(1,-1), columns=sd_stats[i].columns.values)
        sdmax_frame = pd.DataFrame(sd_stats[i].loc["max"].values.reshape(1,-1), columns=sd_stats[i].columns.values)
    else:
        mean_frame = pd.concat((mean_frame, pd.DataFrame(mean_stats[i].loc["mean"].values.reshape(1,-1), columns=mean_stats[i].columns.values)),axis=0)
        meanstd_frame = pd.concat((meanstd_frame, pd.DataFrame(mean_stats[i].loc["std"].values.reshape(1,-1), columns=mean_stats[i].columns.values)),axis=0)
        sd_frame = pd.concat((sd_frame, pd.DataFrame(sd_stats[i].loc["mean"].values.reshape(1,-1), columns=sd_stats[i].columns.values)),axis=0)
        sdmax_frame = pd.concat((sdmax_frame, pd.DataFrame(sd_stats[i].loc["max"].values.reshape(1,-1), columns=sd_stats[i].columns.values)),axis=0)
        
mean_frame["Load"] = x_plot
mean_frame.plot(x = "Load", xlabel=("Load (kN)"), ylabel=("Mean Mean"))
sd_frame["Load"] = x_plot
sd_frame.plot(x = "Load", xlabel=("Load (kN)"), ylabel=("Standard Deviation Mean"))
sdmax_frame["Load"] = x_plot
sdmax_frame.plot(x = "Load", xlabel=("Load (kN)"), ylabel=("Maximum Standard Deviation"))
# below code not so useful - standard deviation of mean depends on load, but only because of component due variation across displacement field
# meanstd_frame["Load"] = x_plot
# meanstd_frame.plot(x = "Load", xlabel=("Load (kN)"), ylabel=("Mean Standard Deviation"))

# mean_mean = [frame.loc["mean"] for frame in mean_stats]
# print(mean_mean)

plt.show()