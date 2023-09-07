from FileSeries import *

def subtract_disp(Files, column_labels = [["x", "u"],["y", "v"],["z","w"]]):
    # Subtract the displacement from the coordinates of each dataset to retrieve the undeformed shape
    # Files = FileSeries class containing info on the location of input csvs and desired output location
    # Define list of coordinate indices, and their corresponding displacement index
    # Column inds = list of coordinate labels and the corresponding displacement component
    Files.read_data()
    # Loop through all files in the folder
    print("Subtracting displacement")
    for i, File in enumerate(Files.files):
        print(i)
        # Subtract each of the displacement components from their corresponding coordinate, and 
        # create a new entry with the result
        for labels in column_labels:
            # Insert new entry to the right of the current column, subtracting the displacement component 
            # from the coordinate
            File.data.insert(loc = File.data.columns.get_loc(labels[0])+1, 
                             column = labels[0]+"_0", 
                             value = (File.data[labels[0]] - File.data[labels[1]]))

    Files.dump()

if __name__ == "__main__":
    folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Alvium Pair 03\\Export_2"
    # folder = "..\\Failure\\Processed DIC Data\\Individual Fields of View\\Manta Camera Pair\\Export_2"
    Files = FileSeries(folder=folder,in_sub_folder="Extracted_qoi",out_sub_folder="Data_subtracted_disp")
    subtract_disp(Files)