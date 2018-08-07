import numpy as np
import pandas as pd
import parmed.unit as u
import os
import sys
import pdb
from multiprocessing import Pool
import itertools

def main():
    all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing)]
    curr_dir = os.getcwd()
    with Pool(8) as p:
        p.starmap(do_stuff, zip(itertools.repeat(curr_dir),all_sweeps))

def do_stuff(curr_dir,sweep_dir):
    os.chdir(os.path.join(curr_dir, sweep_dir))
    sims = [thing for thing in os.listdir() if os.path.isdir(thing) and 'Sim' in thing]
    for sim in sims:
        os.chdir(os.path.join(curr_dir, sweep_dir, sim))

        print((sweep_dir, sim))
        files = [thing for thing in os.listdir() if os.path.isfile(thing) and 'forceout' in thing and 'condensed' not in thing]
        for thing in files:
            condense(thing)

def condense(filename):
    data = np.loadtxt(filename)
    if data.shape[0] > 1:
        timestep =  1 * u.femtoseconds
        timelength = 1*u.nanoseconds
        # First column is the timestep
        last_timestep = 1e6
        condensed_timeseries = []
        condensed_forceseries = []
        is_condensed = False
        row = 0
        while not is_condensed:
            if row == 0:
                condensed_timeseries.append(data[0,0])
                condensed_forceseries.append(data[0,2])
            elif row < data.shape[0]:
                if data[row,0] < last_timestep:
                    if abs(data[row,0] - condensed_timeseries[-1]) > 1:
                        condensed_timeseries.append(data[row,0])
                        condensed_forceseries.append(data[row,2])
                else:
                    is_condensed=True
            else:
                is_condensed =True
            row +=1
        np.savetxt('condensed_{}'.format(filename), np.column_stack((condensed_timeseries,
                                                                condensed_forceseries)))
                
if __name__ == "__main__":
    main()
