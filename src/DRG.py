from utils import *
from cantera import Solution, suppress_thermo_warnings

sys.path.append("dep/pymars")
from argparse import ArgumentParser
from os.path import splitext
from warnings import warn

# local imports
import soln2cti
from drg import run_drg
from sensitivity_analysis import run_sa

# ===========================
# setting
output_name = "DMEsk47"
model_file = "mech/DME/DME2000.cti"
conditions = "mech/DME/conditions.txt"
target_species = ['CH3OCH3','O2']
retained_species = ['CH3OCH3','O2','N2','CO2','H2O']
error = 5.
method = "DRG"
run_sa = False
epsilon_star = 0.01
# thresh = np.array(range(1,100))*0.005
thresh = [0.195]

# ===========================
# preparing
solution_object = Solution(model_file)
final_error = [0]

# running
reduced_model = run_drg(solution_object, conditions, error, target_species,
                        retained_species, model_file, final_error, thresh)
if run_sa:
    print("Running sensitivity analysis")
    reduced_model = run_sa(solution_object, reduced_model, epsilon_star, final_error,
                           conditions, target_species, retained_species, error)

ori_s = [s.name for s in solution_object.species()]
new_s = [s.name for s in reduced_model.species()]
del_s = [s for s in ori_s if s not in new_s]
cprint("\nFinal counts: sp:%d, reac:%d"%(len(reduced_model.species()),len(reduced_model.reactions())),'g')
cprint("Species deleted:"," ".join(del_s),'g')

reduced_model.name = output_name
output_file = soln2cti.write(reduced_model)
