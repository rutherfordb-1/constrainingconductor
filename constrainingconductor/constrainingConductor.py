import os
import pdb
import warnings
import random
import sys
import subprocess
import numpy as np
import mdtraj


class constrainingConductor():
    """ Broad class that orchestrates and executes 
    all stages of Z constraining"""
    def __init__(self, grofile, topfile, z0=1.0, dz=0.2, n_windows=40, n_tracers=8, 
            z_windows=None, auto_detect=True, center=True,
            GMX_CMD="gmx", MDRUN_CMD="gmx mdrun",MPIRUN_CMD="", 
            LMP_CMD="lmp", baseDir=None):

        """
        Parameters
        ----------
        z0 : float
            Initial z coordinate
        dz : float
            window spacing
        n_windows : int
            Number of z-windows
        N_tracer : int 
            Number of tracer molecules
        Z_windows : str
            Filename of z-windows if specified
        auto_detect : Boolean
            If true, automatically generate z windows based on bilayer CoM
        grofile : string
            Filename of gmx structure file
        topfile : string
            Filename of gmx topology file

            """

        self._n_tracers = n_tracers
        self._n_windows = n_windows
        self._grofile = os.path.join(baseDir,grofile)
        self._topfile = os.path.join(baseDir,topfile)
        self._GMX_CMD = GMX_CMD
        self._MDRUN_CMD = MDRUN_CMD
        self._MPIRUN_CMD = MPIRUN_CMD
        self._LMP_CMD = LMP_CMD
        self._dz = dz
        self._z0 = z0
        self._traj = mdtraj.load(self._grofile)
        self._tracers = None
        self._windows = None
        self._n_sims = int(self._n_windows/self._n_tracers)
        self._baseDir = baseDir
        if center:
            self._centerBilayer()

        if auto_detect:
            self._generateWindows()

        elif z_windows is None:
            print("Generating windows from window specifications")
       
            self._windows = list()
            for i in range(self._n_windows):
                   self._windows.append(self._z0 + (i * self._dz)) 
        elif z_windows:
            self._setWindowsFromFile(z_windows)
        else:
            sys.exit("Specify constrainingConductor initializaion parameters!")
        self.selectTracers()

    @property
    def baseDir(self):
        return self._baseDir
    @property
    def grofile(self):
        return self._grofile

    @grofile.setter
    def grofile(self, grofile):
        self._grofile = grofile
        self._traj = mdtraj.load(grofile)

    @property
    def n_sims(self):
        return self._n_sims

    @property
    def n_tracers(self):
        return self._n_tracers

    @property
    def tracers(self):
        return self._tracers

    @tracers.setter
    def tracers(self, tracerfile):
        self._tracers = np.loadtxt(tracerfile, dtype='int')

    @property
    def windows(self):
        return self._windows

    @windows.setter
    def windows(self, windowfile):
        self._setWindowsFromFile(windowfile)


    def writeNdxFile(self, filename=None, pull_group_names=None):
        print("*"*20)
        print("Writing index file...")
        print("*"*20)
        if filename:
            p = subprocess.Popen('echo q | {0} make_ndx -f {1} -o {2}.ndx'.format(self._GMX_CMD, self._grofile, filename), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
        else:
            warnings.warn("Filename unspecified, using gmx default index.ndx", UserWarning)
            p = subprocess.Popen('echo q | {0} make_ndx -f {1}'.format(self._GMX_CMD, self._grofile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()


        if pull_group_names:
            with open(filename+".ndx", "a") as f:
                tracer_atoms = []
                for i, tracer in enumerate(self._tracers):
                    tracer_atoms.append(self._traj.topology.select('resid {}'.format(tracer-1)))

                ##Loop through tracer names and oxygen indices
                for i, groupname in enumerate(pull_group_names):
                    f.write('[ {} ]\n'.format(groupname))
                    f.write('{:10.0f}\t{:10.0f}\t{:10.0f}\n'.format(
                        float(tracer_atoms[i][0]+1), float(tracer_atoms[i][1]+1), float(tracer_atoms[i][2]+1)))



    def _generateWindows(self):
        """ Build simWindows out from the center of mass"""
        print("*"*20)
        print("Generating windows from center of mass...")
        print("*"*20)
        non_water = self._traj.topology.select('not water')
        sub_traj = self._traj.atom_slice(non_water)
        com = mdtraj.compute_center_of_mass(sub_traj)
        center_z = round(com[0,2], 2)
        self._windows = np.zeros(self._n_windows)
        midpoint = int(self._n_windows/2)
        self._windows[midpoint] = center_z

        # Fill in lower half from centerpoint outward
        for i in range(1, midpoint+1):
            self._windows[midpoint-i] = center_z - i * self._dz
        # Fill in top half from centerpoint outward
        for i in range(midpoint+1, self._n_windows):
            self._windows[i] = self._windows[0] + i * self._dz

    def writeWindows(self,filename):
        np.savetxt(filename, self._windows)

    def writeTracers(self,filename):
        np.savetxt(filename, self._tracers, fmt='%i')

    def selectTracers(self):
        """Randomly select tracers"""
        waters = self._traj.topology.select('water and name O')
        atoms = [a for a in self._traj.topology.atoms]
        top_waters = [i for i in waters if self._traj.xyz[0,i,2] <= 
                self._traj.unitcell_lengths[0,2]/3]
        tracer_list = np.random.choice(top_waters, self._n_tracers, replace=False)
        self._tracers =  [atoms[t].residue.index+1 for t in tracer_list]

    def writePullingMdp(self, simWindows, filename="pulling", 
        pull_coord_k=1000, pull_coord_rate=0,
        dt=0.002, t_pulling=1e3, pull_nstfout=5000, pull_nstxout=5000,
        moving=False, z_windows=None,
        T=305, P=1.0): 
        """ Write out MDP parameters for GMX pulling"""



        # Pulling parameters
        pull = 'yes'
        pull_ncoords = self._n_tracers
        pull_ngroups = self._n_tracers
        pull_group_names = ['Tracer'+str(i) for i in self._tracers]
        pull_coord_groups = ['0 {}'.format(i+1) for i in range(self._n_tracers)]
        self.writeNdxFile(filename=filename, pull_group_names=pull_group_names)

        # Pulling origins
        pull_coord_origins = []
        for i, tracer in enumerate(self._tracers):
            xyz = self._getTracerCoordinates(tracer)
            if not moving:
                xyz[2] = simWindows[i]
            pull_coord_origins.append(xyz)

        # Pulling rates
        pull_coord_rate_list = np.zeros(self._n_tracers)
        if moving:
            for i, (tracer, window) in enumerate(zip(self._tracers, simWindows)):
                pull_coord_rate_list[i] = self._calcPullingRate(tracer, window, t_pulling)


        pull_coord_type = 'umbrella'
        if P:
            pull_coord_geometry = 'direction'
        elif not P:
            pull_coord_geometry = 'direction-periodic'
        pull_coord_dim = 'N N Y'
        pull_coord_vec = '0 0 1'
        pull_coord_start = 'no'

        # General MD things
        integrator = 'md'
        nsteps = int(t_pulling/dt)
        comm_mode = 'Linear' 
        nstcomm = 1 
        comm_grps = 'non-water water'
        nstxout = 0 
        nstvout = 0 
        nstfout = 0
        nstxtcout = int(10/dt) 
        nstenergy = int(10/dt) 
        nstlog = int(10/dt) 
        nstcalcenergy = 1


        #Bond parameters 
        constraint_algorithm = 'lincs'
        constraints = 'all-bonds'
        lincs_iter = 1
        lincs_order = 4
        
        #Neighbor searching
        cutoff_scheme = 'Verlet'
        nstlist = 10
        rcoulomb = 1.4
        rvdw = 1.4
        
        #Electrostatics
        coulombtype = 'PME'
        fourierspacing = 0.16
        pme_order = 4
        
        #Temperature coupling
        tcoupl = 'nose-hoover'
        tc_grps = '{:8s}\t{:8s}'.format('non-water', 'water')
        tau_t = '{:8s}\t{:8s}'.format('0.4', '0.4')
        ref_t = '{:8s}\t{:8s}'.format(str(T), str(T))

        #Pressure coupling
        if P:
            pcoupl = 'Parrinello-Rahman'
            pcoupltype = 'semiisotropic'
            tau_p = 2.0
            ref_p = '{} {}'.format(P,P)
            compressibility = '4.5e-5 4.5e-5'
            refcoord_scaling = 'com'
        elif not P:
            pcoupl = 'no'
            pcoupltype = ''
            tau_p = 0.0
            ref_p = '0 0'
            compressibility = '0 0'
            refcoord_scaling = 'com'

        #Misc stuff
        continuation = 'yes'
        gen_vel = 'no'
        pbc = 'xyz'
        DispCorr = 'EnerPres'

        #Writing
        with open(filename+".mdp", 'w') as f:
            f.write('; Run MDP parameters\n')

            f.write('{:25s} = {}\n'.format('integrator',integrator))
            f.write('{:25s} = {}\n'.format('dt', str(dt)))
            f.write('{:25s} = {}\n'.format('nsteps', str(nsteps)))
            f.write('{:25s} = {}\n'.format('comm-mode', str(comm_mode)))
            f.write('{:25s} = {}\n'.format('nstcomm', str(nstcomm)))
            f.write('{:25s} = {}\n'.format('comm-grps', str(comm_grps)))
            f.write('\n; Output parameters\n')
            f.write('{:25s} = {}\n'.format('nstxout', str(nstxout)))
            f.write('{:25s} = {}\n'.format('nstvout', str(nstvout)))
            f.write('{:25s} = {}\n'.format('nstxtcout', str(nstxtcout)))
            f.write('{:25s} = {}\n'.format('nstenergy', str(nstenergy)))
            f.write('{:25s} = {}\n'.format('nstlog', str(nstlog)))
            f.write('{:25s} = {}\n'.format('nstfout', str(nstfout)))
            f.write('{:25s} = {}\n'.format('nstcalcenergy', str(nstcalcenergy)))
            f.write('\n; Bond parameters\n')
            f.write('{:25s} = {}\n'.format('continuation', str(continuation)))
            f.write('{:25s} = {}\n'.format('constraint-algorithm', str(constraint_algorithm)))
            f.write('{:25s} = {}\n'.format('constraints', str(constraints)))
            f.write('{:25s} = {}\n'.format('lincs-iter', str(lincs_iter)))
            f.write('{:25s} = {}\n'.format('lincs-order', str(lincs_order)))
            f.write('\n; Neighbor searching\n') 
            f.write('{:25s} = {}\n'.format('cutoff-scheme', str(cutoff_scheme)))
            f.write('{:25s} = {}\n'.format('nstlist', str(nstlist)))
            f.write('{:25s} = {}\n'.format('rcoulomb', str(rcoulomb)))
            f.write('{:25s} = {}\n'.format('rvdw', str(rvdw)))
            f.write('\n; Electrostatics\n')
            f.write('{:25s} = {}\n'.format('coulombtype', str(coulombtype)))
            f.write('{:25s} = {}\n'.format('fourierspacing', str(fourierspacing)))
            f.write('{:25s} = {}\n'.format('pme_order', str(pme_order)))
            f.write('\n; Temperature coupling\n')
            f.write('{:25s} = {}\n'.format('tcoupl', str(tcoupl)))
            f.write('{:25s} = {}\n'.format('tc_grps', str(tc_grps)))
            f.write('{:25s} = {}\n'.format('tau_t', str(tau_t)))
            f.write('{:25s} = {}\n'.format('ref_t', str(ref_t)))
            f.write('\n; Pressure coupling\n')
            f.write('{:25s} = {}\n'.format('pcoupl', str(pcoupl)))
            f.write('{:25s} = {}\n'.format('pcoupltype', str(pcoupltype)))
            f.write('{:25s} = {}\n'.format('tau_p', str(tau_p)))
            f.write('{:25s} = {}\n'.format('ref_p', str(ref_p)))
            f.write('{:25s} = {}\n'.format('compressibility', str(compressibility)))
            f.write('{:25s} = {}\n'.format('refcoord_scaling', str(refcoord_scaling)))
            f.write('\n; Misc stuff\n')
            f.write('{:25s} = {}\n'.format('gen_vel', str(gen_vel)))
            f.write('{:25s} = {}\n'.format('pbc', str(pbc)))
            f.write('{:25s} = {}\n'.format('DispCorr', str(DispCorr)))
            f.write('\n; Pull parameters\n')
            f.write('{:25s} = {}\n'.format('pull', str(pull)))
            f.write('{:25s} = {}\n'.format('pull-nstxout', str(pull_nstxout)))
            f.write('{:25s} = {}\n'.format('pull-nstfout', str(pull_nstfout)))
            f.write('{:25s} = {}\n'.format('pull-ngroups', str(pull_ngroups)))
            f.write('{:25s} = {}\n'.format('pull-ncoords', str(pull_ncoords)))
            for i in range(self._n_tracers):
                f.write('{:25s} = {}\n'.format('pull-group'+str(i+1)+'-name', pull_group_names[i]))        
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-groups', pull_coord_groups[i]))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-type', pull_coord_type))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-geometry', pull_coord_geometry))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-vec', pull_coord_vec))
                f.write('{:25s} = {:<8.3f} {:<8.3f} {:<8.3f}\n'.format('pull-coord'+str(i+1)+'-origin', 
                    pull_coord_origins[i][0], pull_coord_origins[i][1], pull_coord_origins[i][2]))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-dim', pull_coord_dim))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-rate', pull_coord_rate_list[i]))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-k', pull_coord_k))
                f.write('{:25s} = {}\n'.format('pull-coord'+str(i+1)+'-start', pull_coord_start))




    def _getTracerCoordinates(self, tracer):
        """ 

        Parameters
        ----------
        tracer : int
            Molecule number

        Notes
        -----
        Gromacs is 1-indexed, while mdtraj is 0-indexed
        Getting the coordinates of gmx resid 3 means
        getting the coordinates of mdtraj resid 2
            """
        tracer_atoms = self._traj.topology.select('resid {}'.format(tracer-1))
        sub_traj = self._traj.atom_slice(tracer_atoms)
        com = mdtraj.compute_center_of_mass(sub_traj)

        return [com[0,0], com[0,1], com[0,2]]

    def _calcPullingRate(self, tracer, window, t_pulling):
        xyz = self._getTracerCoordinates(tracer)

        distance_to_traverse = float(window) - xyz[2] #nm
        pull_rate = distance_to_traverse/t_pulling #nm/ps
        return pull_rate

    #def grompp(self, pullingMdp, index, output):
    def grompp(self, filename):
        print("*"*20)
        print("Grompping {}...".format(filename))
        with open('grompp_{}.log'.format(filename), 'w') as f:
            p = subprocess.Popen('{0} grompp -f {1}.mdp -c {2} -p {3} -n {4}.ndx -o {5} -maxwarn 2'.format(
                self._GMX_CMD, filename, self._grofile, 
                self._topfile, filename, filename),
                shell=True, stdout=f, stderr=f)
            p.wait()
        if not os.path.exists(filename+".tpr"):
            print("ERROR: {} not found".format(filename+".tpr"))

        print("*"*20)

    def mdrun(self, output):
        print("*"*20)
        print("MDRunning {}...".format(output))

        with open('mdrun_{}.log'.format(output), 'w') as f:
            p = subprocess.Popen('{0} {1} -deffnm {2}'.format(
                self._MPIRUN_CMD, 
                self._MDRUN_CMD, output), shell=True, stdout=f, stderr=f)
            p.wait()

        if not os.path.exists(output+".gro"):
            print("ERROR: {} not found".format(output+".gro"))

        print("*"*20)

    def lmprun(self, output, lmp_input):
        print("*"*20)
        print("LmpRunning {}...".format(output))

        with open('lmprun_{}.log'.format(output), 'w') as f:
            p = subprocess.Popen('{0} {1} < {2} >& {3}'.format(
                self._MPIRUN_CMD,
                self._LMP_CMD, lmp_input, f.name), shell=True, stdout=f, stderr=f)
            p.wait()
        if not os.path.exists('restartfile'):
            print("ERROR: Lmps simulation failed")
        print("*"*20)

    def _centerBilayer(self):
        print("*"*20)
        print("Centering bilayer...")
        print("*"*20)
        non_water = self._traj.topology.select('not water')
        sub_traj = self._traj.atom_slice(non_water)
        # Get center of mass of the bilayer 
        com = mdtraj.compute_center_of_mass(sub_traj)
        box = [self._traj.unitcell_lengths[0,0], 
                self._traj.unitcell_lengths[0,1], 
                self._traj.unitcell_lengths[0,2]]

        # Shift all coordinates so bilayer is at center of box
        for i, val in enumerate(self._traj.xyz):
            self._traj.xyz[i] = [[x, y, z-com[0,2]+(box[2]/2)] for x,y,z 
                    in self._traj.xyz[i][:]]

        # Go back and fix for pbc
        for i, val in enumerate(self._traj.xyz):
            for j,atom in enumerate(self._traj.xyz[i]):
                for k, coord in enumerate(self._traj.xyz[i,j]):
                    if coord < 0:
                        self._traj.xyz[i,j,k] += box[k]
                    elif coord > box[k]:
                        self._traj.xyz[i,j,k] -= box[k]
        self._traj.save(os.path.join(self._baseDir,'centered.gro'))
        #self.grofile = 'centered.gro'
        self._grofile = (os.path.join(self._baseDir,'centered.gro'))

    def _setWindowsFromFile(self, z_windows):
        print("*"*20)
        print("Loading windows from file...")
        print("*"*20)
        self._windows = np.loadtxt(z_windows)
        self._z0 = round(self._windows[0], 2)
        self._dz = round(self._windows[1] - self._windows[0],2)
        self._n_windows = len(self._windows)
    

