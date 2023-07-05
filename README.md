Contains code used to process point clouds from processed Digital Image Correlation (DIC) data, and align these with Abaqus Finite Element models.
Includes code for:
1) Loading DIC from csv, extracting quantities of interest, and if necessary removing outliers.
2) Rotating/translating DIC coordinates and displacements according to specified translation/rotation matrices for alignment with the model
3) Projecting each point onto the model surface, matching to an element, and determining the natural coordinates within that element
4) Output of processed data to csvs