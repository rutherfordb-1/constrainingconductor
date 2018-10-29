import groToLmps
import pdb
import subprocess
import os
import numpy as np
import mdtraj
import argparse

Lmp_FF_file = "/raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/"

def _write_input_header(f, temp=305.0, Nrun=1000000, Nprint=1000, 
        structure_file='Stage4_Eq0.lammpsdata', 
        tracers=None, tracer_atom_indices=None):
    f.write("""clear
variable Nprint equal {Nprint} 
variable Nrun equal {Nrun} 
variable temperature equal {temp}

units real
atom_style full

pair_style lj/charmm/coul/long 10.0 12.0
bond_style harmonic
angle_style charmm
dihedral_style charmm
improper_style harmonic
special_bonds charmm
kspace_style pppm 1.0e-4
neighbor 2.0 bin 

read_data {structure_file}
write_restart sim.restart.0

clear
read_restart sim.restart.*
kspace_style pppm 1.0e-4
neighbor 2.0 bin 

variable Nprint equal {Nprint} 
variable Nrun equal {Nrun} 
variable temperature equal {temp}

restart 50000 sim.restart.*


    """.format(**locals()))

    water1_type, water2_type, water_bond, water_angle = _parse_water_info(structure_file)


    f.write("""
group water type {0} {1}
group tracers id {2}
group allbuttracers subtract water tracers
group bilayer subtract all water
fix 11 all shake 0.0001 10 10000 b {3} a {4}
""".format(water1_type, water2_type, 
            ' '.join(np.asarray(tracer_atom_indices, dtype=str).flatten()),
            water_bond, water_angle))

def _parse_water_info(structure_file):
    """ Get information about the water molecules
    Look through the lammps data file for the water's two atomtypes, bond type,
    and angle type. The idea is to find a directive and look at the previous 3 or
    4 lines to get the atom/bond/angle type

    Parameters
    ---------
    structure_file : str
        filename of the lammps structurefile

    Returns
    -------
    water1_type : str
    water2_type : str
    water_bond : str
    water_angle : str
        """
    lmp_lines = open(structure_file, 'r').readlines()
    # Find Bonds directive for getting the atomtypes
    bond_line = [line_index for line_index, line in enumerate(lmp_lines) if 'Bonds' in line][0]
    water1_type = lmp_lines[bond_line - 4].split()[2]
    water2_type = lmp_lines[bond_line - 3].split()[2]

    # Find Angles directive for getting bondtype
    angle_line = [line_index for line_index, line in enumerate(lmp_lines) if 'Angles' in line][0]
    water_bond = lmp_lines[angle_line - 2].split()[1]

    # Find Dihedrals directive for getting angletype
    dihedral_line = [line_index for line_index, line in enumerate(lmp_lines) if 'Dihedrals' in line][0]
    water_angle = lmp_lines[dihedral_line - 2].split()[1]

    return water1_type, water2_type, water_bond, water_angle


def _write_rest(f, z_windows, force_indices, tracer_atom_indices,
                record_force=True):
    f.write("""
variable ke equal ke
variable enthalpy equal enthalpy
variable pe equal pe
variable step equal step
variable temp equal temp
variable press equal press
variable vol equal vol
variable time equal time

variable lx equal lx
variable ly equal ly
variable lz equal lz

timestep 1.0

fix 3 all print ${Nprint} "${step} ${time} ${pe} ${press} ${temp} ${lx} ${ly} ${lz}" file system.log screen no
fix 4 all npt temp ${temperature} ${temperature} 10.0 aniso 1.0 1.0 100.0
fix 5 all momentum 1 linear 1 1 1
thermo ${Nprint}
dump d1 all dcd 5000 trajectory.dcd
""")

    for i, (window, force_index, atom_indices) in enumerate(zip(z_windows,
                            force_indices, tracer_atom_indices)):
        f.write("group t{0} id {1}\n".format(i, ' '.join(np.asarray(atom_indices, 
                                                        dtype=str))))

        if record_force:
            f.write("compute tracerfz{0} t{0} reduce sum fz\n".format(i))
            f.write("variable redforce{0} equal c_tracerfz{0}\n".format(i))
            f.write("fix pt{0} all print 3 \"${{step}} ${{time}} ${{redforce{0}}}\" append forceout{1}.dat screen no\n".format(i, force_index))

        f.write("fix 6{0} t{0} recenter NULL NULL {1} shift t{0}\n".format(i, window))
        f.write("fix 7{0} t{0} momentum 1 linear 0 0 1\n".format(i))
        f.write("\n")
    f.write("""
run ${Nrun} upto

write_restart restartfile
""")
   
def _prepare_lmps(eq_structure, z_windows, tracers,
        force_indices, forcefield_files=['foyer_charmm.xml']):
    """ Convert structure to lammps and generate input file"""
    
    traj = mdtraj.load(eq_structure)
    groToLmps.convert(eq_structure, forcefield_files=forcefield_files)
    midplane = np.mean(traj.unitcell_lengths[:,2])
    tracer_atom_indices = _get_atom_indices(traj, tracers)

    ## Load in the gmx z windows, scale/shift appropriately
    #z_windows = [np.round(10*(z - traj.unitcell_lengths[0][2]/2),2) for z in z_windows]
    # Scale and shift gromacs windows accordingly
    # Units are now angstroms and the box still is bottom left origin
    # Correct for PBCs
    new_z_windows = np.zeros_like(z_windows)
    for i, (val, atom_index) in enumerate(zip(z_windows, tracer_atom_indices[::3])):
        new_z_windows[i] = val*10 
        # If the windows are too close to the periodic box boundaries, 
        # the lammps simulation will crash
        # Furthermore, ensure the windows are on the same 'side' as the tracers
        # themselves
        if val < 0.2 or val >= traj.unitcell_lengths[0][2] - 0.2:
            if traj.xyz[0, atom_index-1, 2] < midplane:
                new_z_windows[i] = 2
            else:
                new_z_windows[i] = 10*traj.unitcell_lengths[0][2] - 2
    np.savetxt('z_windows_lmps.out', new_z_windows)
    n_windows = np.shape(new_z_windows)[0]

    ## Load in the tracer information
    n_tracers = np.shape(tracers)[0]



    with open('Stage5_ZCon.input','w') as f:
        _write_input_header(f, structure_file=eq_structure[:-4]+'.lammpsdata', 
                            tracers=tracers, tracer_atom_indices=tracer_atom_indices)
        _write_rest(f, new_z_windows, force_indices, tracer_atom_indices)

def _get_atom_indices(traj, tracers):
    """ Get atom indices correspond to each tracer molecule

    Parameters
    ---------
    traj : mdTraj.Trajectory
    tracers : n x 1 tuple
        List of molecule IDs or residue IDs

    Returns
    ------
    tracer_indices: n x 3 tuple
        List where each element corresponds to each molecule
        Each 3-tuple corresponds to the atom indices in that molecule

    Notes
    -----
    MDtraj does everything 0-indexed, but we want everything 1-indexed
    """
    tracer_indices = []
    for tracer in tracers:
        atoms = [a.index + 1 for a in traj.topology.residue(tracer-1).atoms]
        tracer_indices.append(atoms)

    return tracer_indices




def lmps_conversion(eq_structure, windows, tracers, force_indices):
    """ Main method to do lmps structure conversions and input writing"""
    print("*"*20)
    print("Converting in {}".format(os.getcwd()))
    print("*"*20)
    p = subprocess.Popen("cp {}/*.xml . ".format(Lmp_FF_file), shell=True, 
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    _prepare_lmps(eq_structure, windows, 
                tracers, force_indices)


