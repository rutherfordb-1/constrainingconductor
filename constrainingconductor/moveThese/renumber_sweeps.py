import numpy as np
import subprocess
import os

""" To renumber all the sweep folders
I need to move the old sweep folder to a temp folder
then move the temp folders to new sweep folder names
in order to prevent an overwrite clashes"""
all_temp_folders = []
all_new_folders = []
all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and 'sweep' in thing[0:5]]

for i, old_folder in enumerate(all_sweeps):
    temp_folder = "temp_{}".format(old_folder)
    all_temp_folders.append(temp_folder)

    new_folder = "sweep{}".format(i)
    all_new_folders.append(new_folder)

    p = subprocess.Popen("mv {0} {1}".format(old_folder,temp_folder), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    #print("moving {} to {}".format(old_folder, new_folder))

for temp_folder, new_folder in zip(all_temp_folders, all_new_folders):
    p = subprocess.Popen("mv {0} {1}".format(temp_folder, new_folder), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
