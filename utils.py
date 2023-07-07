# Useful functions called by multiple elements
import os

def make_out_folder(folder_paths):
    # Checks if a directory for writing output already exists, and if not, creates it
    if not(os.path.isdir(folder_paths)):
        os.mkdir(folder_paths)