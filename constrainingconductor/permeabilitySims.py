from constrainingconductor import constrainingConductor
import os
import pdb
import numpy as np
import subprocess
dt = 0.002
sweepStart=0
n_sweeps=30

baseDir = os.getcwd()
GMX_CMD = 'gmx'
MDRUN_CMD = 'mdrun -ntomp 8 -gpu_id 01'

grofile='md_pureDSPC.gro'
topfile='pureDSPC.top'

master = constrainingConductor(grofile,topfile, auto_detect=True,center=True,baseDir=baseDir, GMX_CMD=GMX_CMD, MDRUN_CMD=MDRUN_CMD)
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

            #simWindows = master.windows[sim::master.n_sims]
            master.writePullingMdp(**stageInformation[stage])
            master.grompp(stageInformation[stage]['filename'])
            master.mdrun(stageInformation[stage]['filename'])
            master.grofile = stageInformation[stage]['filename']+".gro"
        os.chdir(baseDir)

    print('-'*10 + 'Finished sweep{}'.format(sweep) + '-'*10)


