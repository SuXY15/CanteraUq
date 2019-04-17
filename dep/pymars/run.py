from argparse import ArgumentParser
from os.path import splitext
from warnings import warn
from pymars import pymars
from convert_chemkin_file import convert
import numpy as np

model = "../mech/DME/chem.cti"
conditions = "../mech/DME/conditions.txt"
error = 5.
method = "DRG"
targets = ['CH3OCH3','O2']
retained_species = ['CH3OCH3','O2','N2','CO2','H2O']
run_sa = False
epsilon_star = 0.01

thresh_arr = np.array(range(1,10))*0.05
thresh_arr = [0.170]

pymars(model, conditions, error, method, targets, thresh_arr,
        retained_species, run_sa, epsilon_star)