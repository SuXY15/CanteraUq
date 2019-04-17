"""Contains main driver function for pyMARS program."""
from cantera import Solution, suppress_thermo_warnings
# local imports
import soln2cti
from drgep import run_drgep
from drg import run_drg
from pfa import run_pfa
from sensitivity_analysis import run_sa
from convert_chemkin_file import convert

# Avoid long warnings from Cantera about thermodynamic polynomials
suppress_thermo_warnings()

def pymars(model_file, conditions, error, method, target_species, thresh=None,
           retained_species=None, run_sensitivity_analysis=False, epsilon_star=0.1
           ):
    """Driver function for pyMARS to reduce a model.

    Parameters
    ----------
    model_file : str
        Cantera-format model to be reduced (e.g., 'mech.cti').
    conditions : str
        File with list of autoignition initial conditions.
    error : float
        Maximum error % for the reduced model.
    method : {'DRG', 'DRGEP', 'PFA'}
        Skeletal reduction method to use.
    target_species: list of str
        List of target species for reduction.
    retained_species : list of str, optional
        List of non-target species to always retain.
    run_sensitivity_analysis : bool, optional
        Flag to run sensitivity analysis after completing another method.
    epsilon_star : float, optional
        Epsilon^* value used to determine species for sensitivity analysis.

    Returns
    -------
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file

    Examples
    --------
    >>> pymars('gri30.cti', conditions_file, 10.0, 'DRGEP', ['CH4', 'O2'], retained_species=['N2'])

    """

    solution_object = Solution(model_file)
    final_error = [0]
    
    if method == 'DRG':
        reduced_model = run_drg(solution_object, conditions, error, target_species, retained_species, model_file, final_error, thresh)
    elif method == 'PFA':
        reduced_model = run_pfa(solution_object, conditions, error, target_species, retained_species, model_file, final_error)
    elif method == 'DRGEP':
        reduced_model = run_drgep(solution_object, conditions, error, target_species, retained_species, model_file, final_error)

    if run_sensitivity_analysis:
        print("Running sensitivity analysis")
        reduced_model = run_sa(solution_object, reduced_model, epsilon_star, final_error, conditions, target_species, retained_species, error)
    
    print("\nFinal counts: sp:%d, reac:%d"%(len(reduced_model.species()),len(reduced_model.reactions())))
    ori_s = [s.name for s in solution_object.species()]
    new_s = [s.name for s in reduced_model.species()]
    del_s = [s for s in ori_s if s not in new_s]
    print("Species deleted:"," ".join(del_s))
    output_file = soln2cti.write(reduced_model)
