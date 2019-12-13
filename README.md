# Organizing a set of pulling simulations (for use with CHARMM force field)  

Observe that some files are "modules" (`.py` files that are meant to be imported,
not run), some are "convenience scripts" (`.py` files that help you stay organized
and manage simulations), some are "scripts" (`.py` files that you will run)

* `constrainingConductor.py` is a python class that contains basic information about a molecular system (gro and top files). From these, and parameters you can pass in, MDP files are written, grompp'd, and simulated from this constrainingConductor class.

* `permeabilitySims.py` is a script that utilizes the constrainingConductor to conduct a set of pulling simulations. This is the one you will likely be running. 
Note the initial variables for gromacs and lammps binaries, as well as
argparser arguments for detailing the constraining simulations

* `lmpsUtils.py` is a module that contains functions for writing lammpsdata files, 
lammps input scripts. 
Note the path necessary for `foyer_charmm.xml`, 
which can be found in the mccabegroup git repo.
Note the filenaming for dumping out the various force timeseries.

* `moveThese/` is a folder with some scripts you will need to move to conduct your 
own permeability simulations.

* `moveThese/checkFails.py` is merely a convenience script to iterate through the 
permeability sweeps and simulations, checking for presence of a file 
`restartfile` (or whatever you want). 
This can be useful to figure out if you got to a certain step within the simulation or workflow

* `moveThese/condense_timeseries.py` is a convenience script. 
When running these lammps simulations, they can and will often fail. 
Using lammps restart files, the simulations can be continued and 
force timeseries files will be appended to. 
However, let's say a simulation had run steps 1-1001, 
dumped a restart on step 900, and then crashed on step 1001.
The continued-simulation will start from step 900, and dump forces from step 900 
onwards, even though the original simulation had dumped forces form steps 1-1001.
Condensing the time series will attempt to adjust for these "overlapped" timesteps.
Additionally, even though the lammps input files designates writing out 3 columns: 
timestep, time, force-in-z, condensing
will only write time, force-in-z

* `moveThese/groToLmps.py` is a module that facilitates conversion of `.gro` files to
`.lammpsdata` files. 
An mbuild-recipe for a bilayer is used to create an mbuild bilayer, whose chemical
species and molecular information must match the ordering and contents of the `gro`
file. The mbuild bilayer has its coordinates then updated with the specified `gro`
file. Note: we do not merely load the `gro` file into an mbuild compound
because the `gro` file's particle names do not match the atomtypes in 
`foyer_charmm.xml`, and foyer atom-typing will fail. 
Thus, we must construct an mbuild bilayer from mbuild compounds whose particle 
names will actually get atom-typed in foyer (using the underscore naming
convention to facilitate atom-typing). 
As it is currently, `groToLmps` will construct 64 DSPC, 1280 SOL, 64 DSPC, 
and 1280 SOL in that order

* `moveThese/mass_stage5.py` is a convenience script that builds SLURM scripts 
for NERSC Edison to run stage 5 using Lammps. 

* `moveThese/renumber_sweeps.py` is a convenience script that renumbers
sweeps in a sequential order. If you had sweeps 0, 1, 5, 11, 20, they would map 
to 0, 1, 2, 3, 4
