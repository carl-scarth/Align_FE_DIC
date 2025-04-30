# Align_FE_DIC

## Overview 

Aimed at manipulating experimental point cloud data for coordinate transformations during alignment with Finite Element (FE) models, and expressing points in local element coordinate systems. Expressing points in a common coordinate system enables direct comparison of outputs. This repository was developed for comparing Digital Image Correlation (DIC) displacement data against [ABAQUS](https://www.3ds.com/products/simulia/abaqus) model output using 4-noded (S4) elements. 

## Features

Align_FE_DIC is built to apply a series of manipulations to time-dependent point cloud data. The time-dependency of the point cloud is represented by storing data across a set of *.csv* or *.vtk* files, with filenames suffixed by integers, where increasing time is represented by increasing suffix value. Align_FE_DIC is constructed to apply methods across all spatial points and time indices to produce output with consistent data structure. Data is stored as pandas DataFrames within *FileSeries* objects, which adapt the functionality of pandas to suit the specified data structure.

Key features include:
1) Basic data manipulations including extracting data from specific columns, renaming columns, down-sampling in time and space, and filtering using data values and Boolean comparators.
2) Coordinate shifts and rotations.
3) Projecting points onto the FE mesh.
4) Determining equivalent point cloud natural coordinates within quadratic isoparametric elements.
5) Interpolating FE nodal output to point cloud locations.
6) Inverse mapping of point cloud measurements onto nodes of the FE mesh.

## How to use

### Overview

Align_FE_DIC is constructed to allow a series of point cloud operations to be performed via a high-level python script. 

For a quick-start demonstration of key functionality, `manage_workflow.py` is included as an example of how to set up this script, along with an example FE mesh and point cloud data. Inputs and output point clouds are stored under parent directory "input_output", in subdirectories "input" and "output" respectively. Mesh data and other inputs are stored under the parent directory. This example script may be run to test functionality and modified to suit the required application and data structure. 

Individual methods may also be used as standalone procedures by modifying their `if __name__ == "__main__":` block and running directly from the command line.

More specific inforrmation on usage is detailed below.

### Point cloud input/output data format

Point cloud data is stored using the *FileSeries* class. This object is set up to automatically detect all point cloud (*.csv* or *.vtk*) files within a specified input folder, read and track data via a list of Pandas DataFrames, append/insert outputs of point cloud operations to specified columns, and write the processed data to a specified output folder with consistent structure to the inputs. 

The directory structure of the input and outputs is defined by passing the following optional keyword arguments when initialising a *FileSeries* object, all of which must be character strings:
- **folder**: Parent directory of input and output data. If not specified, taken as the current working directory.
- **in_subfolder**: Subdirectory which contains input *.csv* or *.vtk* files. If not specified, taken as the parent directory.
- **out_subfolder**: Subdirectory to which outputs will be written. If not specified, taken as in_subfolder, and input data is overwritten.
- **del_sub_folder**: Subdirectory to which deleted data points (due to filtering or dropping points with NA output) are written. If not specified, these points are not retained.

### Finite Element model input/output format

Mesh data is stored via the *SurfaceMesh* class, which contains node and element definitions, and methods to calculate centroids and normals. This object is initialised either by passing NumPy arrays of nodal coordinates and element connectivities, or string(s) used to identify *.csv* files containing this information, including the directory. If using *.csv* input the first row is assumed to be a header, and skipped. 

The SurfaceMesh object must be initialised by passing one of the following combinations of keyword arguments:
- **nodes**, **connectivity** (*NumPy arrays*): n_nodes x 3 array of (*float*) nodal coordinates and n_elements x 4 (*integer*) array of element connectivities.
- **file_string** (*string*): String used to identify both *.csv* files containing the nodal and element definitions in format *string+"_nodes.csv"*, *string+"_elements.csv"*.
- **node_file**, **el_file** (*strings*): Names of *.csv* files containing nodal and element definitions respectively.

## Methods of FileSeries class

Basic data processing operations may be performed by calling methods of the FileSeries class. The following methods are included:
- `read_data`: Read data from all input files in the specified input folder.
- `down_sam`: Down-sample point cloud data at regular time index.
- `truncate_data`: Truncate point cloud data up to a maximum time index.
- `extract_qoi`: Extract Quantities of Interest from columns of the point cloud data, and rename if necessary.
- `filter_by_cond`: Filter the data to delete points based on a list of logical conditions applied to the DataFrame columns. Each condition is passed as a string comprised of a column label, comparative operator, and value agaist which the data is compared. Boolean operators can be used to create compound logical conditions. Each term must be separated by a space for parsing (with the exception of defining a negative number). See `manage_workflow.py` for examples.
- `update_datatype`: Change data type of specified column.
- `apply_func_to_data`: Apply a function to selected columns of the DataFrame, and insert outputs at a specified location. Used internally by other standalone methods.
- `dump`: Write data to the specified output directory.
- `dump_data_del`: Write all data which is marked for deletion to del_sub_folder.
- `dump_nas`: Write all data with NAs in any column to del_sub_folder.

### Standalone methods

The main point cloud coordinate transformation and FE model comparison functionality is implemented in standalone functions, which can either by imported into the high-level workflow script, or run directly from the command line. The following methods are included:
- `downsample_data`: Down-sample to reduce the number of spatial points in the cloud.
- `subtract_displacement`: Subtract measured displacements from the coordinates of each point. Useful for back-calculating the undeformed DIC geometry.
- `rename_files`: Rename files to remove common suffixes from MatchID, e.g. "*_0.tiff*". Helps paraview identify *.csvs* as a file series.
- `transform_coords`: Apply coordinate transformations (shifts and rotations) using rotation matrix (R) and/or translation vector (T). These inputs can be loaded from a text file using `transmat_from_file`. Text file format matches those outputted by [CloudCompare](https://www.danielgm.net/cc/) (see dependencies).
- `get_nat_coords`: Project an aligned point cloud onto an FE mesh of quadratic elements, determine equivalent element natural coordinates for each point, and a list of matching element indices.
- `project_points`: Project an aligned point cloud onto an FE mesh without finding natural coordinates. (Note: This is less accurate than `get_nat_coords`, where it is called internally).
- `interp_nodes_to_point`: Interpolate FE nodal quantities (e.g displacements, coordinates) to point cloud locations given by element indices and natural coordinates.
- `fit_point_to_node`: Inverse mapping of point cloud quantities (e.g. displacements, coordinates) at given element indices and natural coordinates onto the FE nodes using least-squares.

## Dependencies

The methodolgy implicitly assumes that displacment varies within each element as described by the shape functions of an [ABAQUS](https://www.3ds.com/products/simulia/abaqus) 4-noded (S4) shell element. The methodology is also applicable to other Finite Element software and element types, although direct pointwise equivalence may be inexact due to different interpolation functions.

The method used to align the point cloud to the FE model assumes a known rotation matrix (R) and translation vector (T), however, no method is provided for determining transformations which give optimal alignment. The fine-registriation capability of CloudCompare is recommended for performing this alignment, and finding R and T. See:  
<https://www.danielgm.net/cc/>

Paraview is recommended for visualising outputs, and alignment with the FE mesh:  
<https://www.paraview.org/>

Python scripts were implemented and tested using Python 3.13.1, and:
- pandas 2.2.3
- numpy 2.2.2
- natsort 8.4.0
- vtk 9.4.1

## To do

- Update `interp_nodes_to_cloud`, and `point_to_nodes` to *SurfaceMesh* syntax.
- Full-field output comparison via residuals.
- 8-noded shells elements.
- Track progress via progress bar.