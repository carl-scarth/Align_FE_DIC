# Useful functions called by multiple elements
import os
import vtk
import pandas as pd

def make_out_folder(folder_paths):
    # Checks if a directory for writing output already exists, and if not, creates it
    if not(os.path.isdir(folder_paths)):
        os.mkdir(folder_paths)

def vtk_to_pandas(vtk_file_path):
    # Read the VTK file
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(vtk_file_path)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()

    # Get the VTK data
    vtk_data = reader.GetOutput()
    # print(dir(vtk_data))
    # Extract data from the VTK dataset
    num_points = vtk_data.GetNumberOfPoints()
    num_arrays = vtk_data.GetPointData().GetNumberOfArrays()
    num_fields = vtk_data.GetFieldData().GetNumberOfArrays()

    # Create a dictionary to store the data
    data_dict = {}

    # Extract point coordinates
    points = vtk_data.GetPoints()
    coordinates = [points.GetPoint(i) for i in range(num_points)]
    data_dict['X'] = [point[0] for point in coordinates]
    data_dict['Y'] = [point[1] for point in coordinates]
    data_dict['Z'] = [point[2] for point in coordinates]

    # Extract scalar and vector data
    for i in range(num_arrays):
        array_name = vtk_data.GetPointData().GetArrayName(i)
        array = vtk_data.GetPointData().GetArray(i)
        num_components = array.GetNumberOfComponents()

        # If the array has only one component, it's a scalar
        if num_components == 1:
            data_dict[array_name] = [array.GetValue(j) for j in range(num_points)]
        else:
            # If the array has more than one component, it's a vector
            for component in range(num_components):
                component_name = f"{array_name}_{component}"
                data_dict[component_name] = [array.GetComponent(j, component) for j in range(num_points)]

    # Extract field data
    field_data = vtk_data.GetFieldData()
    for i in range(num_fields):
        field_array = field_data.GetArray(i)
        field_name = field_array.GetName()
        # In this case, field data only contains scalar values
        data_dict[field_name] = field_array.GetValue(0)
        
    # Create a Pandas DataFrame from the dictionary
    df = pd.DataFrame(data_dict)

    return df