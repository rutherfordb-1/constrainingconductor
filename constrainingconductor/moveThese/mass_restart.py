import os
import numpy as np
import subprocess
import pdb

comment_these = [5, 6, 8,9,10,11,12,13,15,17,21,22,23,24,26,39,54,61,68,75,82,89,96]

def _write_sbatch(sbatch_script, sweepnum, simnum):
    with open(sbatch_script, 'w') as f:
        f.write("""#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -J ay-{0}-{1}
#SBATCH -o ay-{0}-{1}.out
#SBATCH -e ay-{0}-{1}.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
#SBATCH -L SCRATCH
module load gromacs/5.1.4-2

module load lammps/2017.03.31
export OMP_NUM_THREADS=1
srun -n 48 -c 1 lmp_edison < restart.input &> restart.out""".format(sweepnum, simnum))


def _modify_script(restart_script, original_script):
    with open(original_script, 'r') as f:
        original_lines = f.readlines()
    restart_lines = original_lines
    for line in comment_these:
        restart_lines[line] = "#" + restart_lines[line]
    restart_lines[0] = "read_restart restartfile\n"
    restart_lines[2] = "variable Nrun equal 620000\n"
    restart_lines[106] = "write_restart restartfile2\n"
    with open(restart_script, 'w') as restart_file:
        restart_file.write("".join(restart_lines))

if __name__=="__main__":
    curr_dir = os.getcwd()
    all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and
            'sweep' in thing[0:5]]
    for sweep in all_sweeps:
        os.chdir(os.path.join(curr_dir, sweep))
        all_sims = [thing for thing in os.listdir() if os.path.isdir(thing) and 'Sim' in thing]
        
        for sim in all_sims:
            os.chdir(os.path.join(curr_dir, sweep, sim))
            forceout_files = [blah for blah in os.listdir() if os.path.isfile(blah) and 'forceout' in blah]

            if os.path.isfile("Stage5_ZCon.input"):
                print(os.getcwd()+" success")
                _modify_script("restart.input", "Stage5_ZCon.input")
                _write_sbatch("restart.sbatch", sweep, sim)
                p=subprocess.Popen("sbatch restart.sbatch", shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p.wait()
                if len(forceout_files)  <= 4:
                    print("PROBLEM")
                    print(os.getcwd())
            else:
                print(os.getcwd()+" failed")
        os.chdir(curr_dir)
