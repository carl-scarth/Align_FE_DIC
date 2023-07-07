import os
from utils import *

class FileSeries:
    # Contains properties of a series of .csv files associated with a particular test
    def __init__(self,folder="", in_sub_folder="", out_sub_folder=""):
        # Perhaps pass folder as a kwargs? If not passed assume that data is in the working directory
        self.folder = folder # Parent folder in which the data is stored
        print(self.folder)
        self.in_sub_folder = in_sub_folder # Subfolder in which the data is stored
        if out_sub_folder == "":
            self.out_sub_folder = self.in_sub_folder
        else:
            self.out_sub_folder = out_sub_folder
        self.get_paths()
        print(self.in_path)
        print(self.out_path)
        # Create a directory for writing output if required
        make_out_folder(self.out_path)

    def get_paths(self):
        # Get full path to which data is contained
        wd = os.path.abspath("") # Get working directory
        self.in_path = os.path.join(wd,self.folder,self.in_sub_folder)
        self.out_path = os.path.join(wd,self.folder,self.out_sub_folder)
