#Organizing a set of pulling simulations

constrainingConductor is a python class that contains basic information about a molecular system (gro and top files). From these, and parameters you can pass in, MDP files are written, grompp'd, and simulated from this constrainingConductor class.

permeabilitySims.py is a script that utilizes the constrainingConductor to conduct a set of pulling simulations
