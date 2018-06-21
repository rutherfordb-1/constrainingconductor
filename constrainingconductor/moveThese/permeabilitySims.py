from constrainingconductor import constrainingConductor
from constrainingconductor import lmpsUtils
import os
import pdb
import numpy as np
import subprocess
import argparse 

dt = 0.002
parser = argparse.ArgumentParser()
parser.add_argument('-n', dest='n_sweeps', action='store', type=int, default=10)
parser.add_argument('-s', dest='sweepStart', action='store', type=int, default=0)
parser.add_argument('-c', dest='grofile', action='store', type=str, default=None)
parser.add_argument('-p', dest='topfile', action='store', type=str, default=None)
args = parser.parse_args()

sweepStart = args.sweepStart
n_sweeps = args.n_sweeps
stage5_lmps = True

baseDir = os.getcwd()
GMX_CMD = 'gmx'
MDRUN_CMD = 'gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01'
MPIRUN_CMD = ''
LMP_CMD = "lmp_accre"

#grofile='md_pureDSPC.gro'
#topfile='pureDSPC.top'
grofile = args.grofile
topfile = args.topfile
#grofile = '11_DSPC_C18OH_Remco.gro'
#topfile = '11_DSPC_C18OH_Remco.top'

master = constrainingConductor(grofile,topfile, auto_detect=True,
        center=True,baseDir=baseDir, GMX_CMD=GMX_CMD, MDRUN_CMD=MDRUN_CMD, 
        MPIRUN_CMD=MPIRUN_CMD, LMP_CMD=LMP_CMD)
originalGrofile = master.grofile
master.writeWindows('z_windows.out')
for sweep in range(sweepStart, sweepStart+n_sweeps):
    print('-'*10 + 'Beginning sweep{}'.format(sweep) + '-'*10)
    master.grofile = originalGrofile
    master.selectTracers()

    for sim in range(master.n_sims):
        master.grofile = originalGrofile
        sim_folder = 'sweep{}/Sim{}'.format(sweep,sim)
        p = subprocess.Popen('mkdir -p {}'.format(sim_folder), shell=True, 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.wait()
        os.chdir(sim_folder)
        master.writeTracers('tracers.out')
        master.writeWindows('z_windows.out')
        master.writeWindows(os.path.join(baseDir, 'sweep{}/z_windows.out'.format(sweep)))


        for stage in range(1,6):
            stageInformation = { 
                    1: {'filename':"Stage1_Weak{}".format(sim), 'pull_coord_k':40, 
                        'simWindows':master.windows[sim]*np.ones(master._n_tracers)}, 
                    2: {'filename':"Stage2_Strong{}".format(sim), 'pull_coord_k':500, 
                        'simWindows':master.windows[sim]*np.ones(master._n_tracers)}, 
                    3: {'filename':"Stage3_Moving{}".format(sim), 
                        'pull_coord_k':500, 'moving':True, 'P':None,
                        'simWindows': master.windows[sim::master.n_sims]
                        }, 
                    4: {'filename':"Stage4_Eq{}".format(sim), 'pull_coord_k':1000, 
                        'simWindows': master.windows[sim::master.n_sims]},
                    5: {'filename':"Stage5_ZCon{}".format(sim), 'pull_coord_k':1000, 
                        'pull_nstfout':int(0.01/dt), 
                        'simWindows': master.windows[sim::master.n_sims]
                        } 
                    }

            if stage==5 and stage5_lmps:
                p = subprocess.Popen("rm forceout*", shell=True, 
                        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                p.wait()
                p = subprocess.Popen("rm trajectory.lammps", shell=True, 
                        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                p.wait()
                p = subprocess.Popen("rm tracerpos.xyz", shell=True, 
                        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                p.wait()

                force_indices = [sim + int(i*len(master.windows)/master.n_tracers) 
                        for i in range(master.n_tracers)]
                lmpsUtils.lmps_conversion(stageInformation[4]['filename']+".gro",
                        master.windows[sim::master.n_sims],
                        master.tracers, force_indices)

                master.lmprun(stageInformation[5]['filename']+'_lmps',
                        'Stage5_ZCon.input')

            else:
                master.writePullingMdp(**stageInformation[stage])
                master.grompp(stageInformation[stage]['filename'])
                master.mdrun(stageInformation[stage]['filename'])
                try:
                    master.grofile = stageInformation[stage]['filename']+".gro"
                except (IOError,OSError) as e:
                    print("{}.gro not found, ending sweep{}/sim{}".format(stageInformation[stage]['filename'],
                        sweep, sim))
                    break
                p = subprocess.Peopen("rm \"#\"*", shell=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p.wait()
        os.chdir(baseDir)

    print('-'*10 + 'Finished sweep{}'.format(sweep) + '-'*10)


