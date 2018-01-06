import os
import sys
import numpy as np
""" Look for Stage5 pull_f.xvg files in every sim in every sweep,
if not present then this simulation failed somewhere"""
sweeps = sorted([filename for filename in os.listdir() if 'sweep' in filename and os.path.isdir(filename)])


# Read in the forces files, splitting 
# Them into different force files
current_dir = os.getcwd()
for sweep in sweeps:
    os.chdir(os.path.join(current_dir, sweep))
    N_sims = len([filename for filename in os.listdir() if 'Sim' in filename and os.path.isdir(filename)])
    for i in range(N_sims):
        #filename = "Stage5_ZCon"+str(i)+"_pullf.xvg"
        filename = "restartfile"
        os.chdir(os.path.join(current_dir, sweep, "Sim{}".format(i)))
        if not os.path.isfile(filename):
            print("Sweep: " + str(sweep)+', '+ "Sim: " + str(i))



