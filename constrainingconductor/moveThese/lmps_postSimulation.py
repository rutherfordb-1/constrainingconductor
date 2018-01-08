import numpy as np
import subprocess
from multiprocessing import Pool
import os
import itertools
import glob

def copy_force_files(curr_dir, sweep):
    """ Thread function to copy files"""
    os.chdir(os.path.join(curr_dir, sweep))
    sims = [sim for sim in os.listdir() if os.path.isdir(sim) and 'Sim' in sim]
    for sim in sims:
        os.chdir(os.path.join(curr_dir, sweep, sim))
        files = glob.glob('forceout*')
        for thing in files:
            data = np.loadtxt(thing, comments=["#", "@"])
            if np.shape(data)[0] < 1:
                print(os.getcwd())
            #if ".dat" in thing:
            #    p = subprocess.Popen('cp {} {}'.format(thing, 
            #        os.path.join(curr_dir, sweep, thing[:-4])), shell=True, 
            #        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #    p.wait()
            #else:
            #    p = subprocess.Popen('cp {} {}'.format(thing, 
            #        os.path.join(curr_dir, sweep, thing)), shell=True, 
            #        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #    p.wait()
            p = subprocess.Popen('cp forceout* {}'.format( 
                os.path.join(curr_dir, sweep)), shell=True, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
            p = subprocess.Popen("rm \"#\"*", shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE)
            p.wait()


    os.chdir(os.path.join(curr_dir, sweep))


curr_dir = os.getcwd()
sweeps = [sweep for sweep in os.listdir() if os.path.isdir(sweep) and 'sweep' in sweep]
with Pool() as pool:
    pool.starmap(copy_force_files, zip(itertools.repeat(curr_dir),sweeps))

