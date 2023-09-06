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
        # Get all filenames in output folder
        self.get_filenames()
        # Possibly also join with in path

    def get_paths(self):
        # Get full path to which data is contained
        wd = os.path.abspath("") # Get working directory
        self.parent_path = os.path.join(wd,self.folder)
        self.in_path = os.path.join(wd,self.folder,self.in_sub_folder)
        self.out_path = os.path.join(wd,self.folder,self.out_sub_folder)

    def get_filenames(self):
        # Get the name of all csv files in the input folder
        self.in_filenames = [file for file in os.listdir(self.in_path) if ".csv" in file]
        # out filenames will be updated if "remove_suffix" is called
        self.out_filenames = self.in_filenames