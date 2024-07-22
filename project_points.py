import numpy as np
import warnings
from SurfaceMesh import *
from FileSeries import *

def project_points(Files, Mesh, coord_labels = ["x","y","z"], in_sub = [], out_sub = "proj", **kwargs):

    warnings.warn('Points will be projected onto the element with the nearest centroid which may not be correct. For more accurate results use "get_nat_coords" instead')
    # Determine input labels  
    if in_sub:
        in_labels = ["_".join((label,in_sub)) for label in coord_labels]
    else:
        in_labels = coord_labels
    # Apply projection to data
    Files.apply_func_to_data(lambda x:loop_projection(x, Mesh, in_labels), coord_labels, in_sub = in_sub, out_sub = out_sub, message = "Projecting points onto mesh", **kwargs)

def loop_projection(cloud_data, Mesh, coord_labels):
    cloud_xyz = cloud_data[coord_labels].to_numpy() # coordinates
    # Initialise outputs
    xyz_proj = np.empty([cloud_xyz.shape[0],3]) # xyz coordinate of projected point
    xyz_proj.fill(np.nan)
    for i, point in enumerate(cloud_xyz):
        # Find element with the closest centroid to the point, to use as an intial guess for the element
        min_el = find_closest_centroid(point, Mesh.centroids)
        min_el = min_el[0] # Take first entry as we aren't searching over the other nearby elements
        # Project onto the selected element
        point_proj = proj_point_on_element(point, Mesh.elements[min_el].centroid, Mesh.elements[min_el].n)
        xyz_proj[i,:] = point_proj
    return(xyz_proj)

def proj_point_on_element(point, cen, n):
    # Project point onto element with centroid "cen" and unit normal "n"
    # Calcuate vector between the point and element centroid
    v = point - cen
    # Project this vector onto the element normal to get the normal distance
    norm_dist = np.sum(v*n)
    # Move this distance along the normal such that the point lies on the surface
    point_proj = point - norm_dist*n
    return(point_proj)

def find_closest_centroid(point, centroids, n_sort = 1):
    # Find the index df element with closest centroid to point "point",
    # where "centroids" is a numpy array where each row is the centroid of an element
    # Find the distance from each point to the centroid of all elements
    # if sorted_points is passed, instead returns the closet n_sort points
    cen_dist = np.linalg.norm(point - centroids, axis=1)
    min_ind = np.argsort(cen_dist)[0:n_sort]
    # Negligable difference in runtime to output multiple minima so just do that by default
    # min_ind = np.argmin(cen_dist)

    return(min_ind)

if __name__ == "__main__":
    # Create file series and load in data
    folder = "E:\\MengYi_Data\\CS02P_DIC\\Right Camera Pair"
    Files = FileSeries(folder=folder,in_sub_folder="Trimmed_Data", out_sub_folder="Projection_Test")

    # Load in mesh and create mesh object, containing nodal coordinates
    # and connectivities, as well as methods for calculating centroids,
    # normals etc
    node_file = "E:\\MengYi_Data\\coords_undeformed.csv"
    el_file = "E:\\MengYi_Data\\element_quad.csv"
    # Construct mesh object based on connectivities, and calculate 
    # element normals and centroids
    Mesh = SurfaceMesh(from_file = True, node_file=node_file, el_file=el_file)

    # Project DIC onto mesh surface
    Files.read_data()
    project_points(Files, Mesh, in_sub = "rot")
    Files.dump()