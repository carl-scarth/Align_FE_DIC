import os
from utils import *

class FileSeries:
    # Contains properties of a series of .csv files associated with a particular test
    def __init__(self,folder="", in_sub_folder="", out_sub_folder=""):
        # Perhaps pass folder as a kwargs? If not passed assume that data is in the working directory
        self.folder = folder # Parent folder in which the data is stored
        self.in_sub_folder = in_sub_folder # Subfolder in which the data is stored
        if out_sub_folder == "":
            self.out_sub_folder = self.in_sub_folder
        else:
            self.out_sub_folder = out_sub_folder

        self.get_paths()
        # Create a directory for writing output if required
        make_out_folder(self.out_path)
        # Create a list of all files in the in_sub_folder
        self.get_files()

    def get_paths(self):
        # Get full path to which data is contained
        wd = os.path.abspath("") # Get working directory
        self.parent_path = os.path.join(wd,self.folder)
        self.in_path = os.path.join(wd,self.folder,self.in_sub_folder)
        self.out_path = os.path.join(wd,self.folder,self.out_sub_folder)

    def get_files(self):
        # Create a list of SingleFile objects
        self.files = [SingleFile(filename, self.in_path, self.out_path) for filename in os.listdir((self.in_path))]

    def dump(self):
        # For just renaming it's quicker to copy. Incorporate this here when ready...
        # Write processed data to csv files
        asds

class SingleFile:
    # Properties of an individual file within a FileSeries
    def __init__(self, filename, in_path, out_path):
        self.in_filename = filename
        self.out_filename = filename # Default output filename is the same as the input
        self.in_path = in_path
        self.out_path = out_path
        self.src = os.path.join(self.in_path, filename)
        self.dst = os.path.join(self.out_path, filename)
        self.data = [] # add data from csv files

    def update_dst(self):
        # Update the destination address if this is renamed
        self.dst = os.path.join(self.out_path, self.out_filename)