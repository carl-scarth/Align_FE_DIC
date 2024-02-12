import vtk
import pandas as pd

# Example usage
vtk_file_path = 'Raw_Data\\frame_300.vtk'
df = vtk_to_pandas(vtk_file_path)

# Display the DataFrame
print(df.head())
# df.to_csv(vtk_file_path.replace("vtk","csv"),sep=",",index=False)