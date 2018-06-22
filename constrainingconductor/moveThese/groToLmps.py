import numpy as np
import pdb
import mdtraj
import mbuild as mb
from mbuild.formats.lammpsdata import write_lammpsdata
from foyer import Forcefield

# Import statements for molecule prototypes
import atomistic.dppc.DPPC as DPPC
import atomistic.dspc.DSPC as DSPC
import atomistic.c12oh.oh12 as oh12
import atomistic.c16oh.oh16 as oh16
import atomistic.c24oh.oh24 as oh24
import atomistic.c12ffa.ffa12 as ffa12
import atomistic.c16ffa.ffa16 as ffa16
import atomistic.c24ffa.ffa24 as ffa24
import atomistic.tip3p_pppm.SOL as SOL


#if __name__ == "__main__":
def convert(correct_structure, 
            forcefield_files=['foyer_charmm.xml', 'foyer_water.xml']):
    """ Convert gmx to lmp structure via foyer

    Parameters
    ---------
    correct_structure : str
        filename of the correct structure we are pulling xyz from
    forcefield_files :
        To be passed to foyer

    Notes
    -----
    We are making a 'fake' mbuild compound and then updating the coordinates. 
    The fake mbuild compound is then converted to parmed structure and then atomtyped
    It is *necessary* to update how we are making this fakae mbuild compound
    It will help to look at the gmx top file to see the order of the molecules
    Remember to specify `use_atom_name` to false in order for the ffxml to succesfully
    ientify atoms in the mb Compound and pmd Structure
        """
    system = mb.Compound()
    for i in range(72):
        system.add(DPPC.DPPC(use_atom_name=False))
    
    for i in range(2160):
        system.add(SOL.SOL(use_atom_name=False))
    system.update_coordinates(correct_structure)
    
    
    # In order to avoid using smarts, define custom elements in parmed
    # by adding underscores to mb particle names
    for num, i in enumerate(system.particles()):
        i.name = "_{}".format(i.name)
    
    structure = system.to_parmed(box=system.boundingbox, residues=set([p.parent.name for p in system.particles()]))
    
    ff = Forcefield(forcefield_files=forcefield_files)
    structure = ff.apply(structure, assert_dihedral_params=False)
    
    # Because mbuild compounds don't pass charges to parmed structures, need to
    # manuallly set the charges AFTER the force field has been applied
    for i, j in zip(system.particles(), structure.atoms):
        j.charge = i.charge
    
    write_lammpsdata(structure, correct_structure[:-4]+'.lammpsdata')
