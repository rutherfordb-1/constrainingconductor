import os
import numpy as np
import subprocess
import pdb

def _write_sbatch(sbatch_script, sweepnum, simnum):
    with open(sbatch_script, 'w') as f:
        f.write("""#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -J c16-{0}-{1}
#SBATCH -o c16-{0}-{1}.out
#SBATCH -e c16-{0}-{1}.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.h.yang@vanderbilt.edu
#SBATCH -L SCRATCH
module load gromacs/5.1.4-2

module load lammps/2017.03.31
export OMP_NUM_THREADS=1
srun -n 48 -c 1 lmp_edison < Stage5_ZCon.input &> Stage5.out""".format(sweepnum, simnum))


if __name__=="__main__":
    curr_dir = os.getcwd()
    all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and
            'sweep' in thing[0:5]]
    for sweep in all_sweeps:
        os.chdir(os.path.join(curr_dir, sweep))
        all_sims = [thing for thing in os.listdir() if os.path.isdir(thing) and 'Sim' in thing]
        
        for sim in all_sims:
            os.chdir(os.path.join(curr_dir, sweep, sim))
            #forceout_files = [blah for blah in os.listdir() if os.path.isfile(blah) and 'forceout' in blah]

            if os.path.isfile("Stage5_ZCon.input"):
                _write_sbatch("stage5.sbatch", sweep, sim)
                p=subprocess.Popen("sbatch stage5.sbatch", shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p.wait()
                print(os.getcwd()+" success")
            else:
                print(os.getcwd()+" failed")
        os.chdir(curr_dir)
