##############
# This function takes the conditions array and converts it into an array of simulation species_objects
#
# conditions array: array of initial conditions read in from the ic file
#  model: model to base the simulation off of
#
# simulation array: an array of unran simulation objects from simulation.py
##############

from simulation import Simulation
import numpy as np
import time

t0 = time.time()
def progress(percent, msg='', show_time=0, width=40):
    '''
    Logging with progress bar
    '''
    percent = 1.0 if percent > 995e-3 else percent+5e-3
    time_str='|%6.1fm|'%((time.time()-t0)/60)
    show_str=('[%%-%ds]' %width) %(int(width * percent)*"#")
    show_str = msg+time_str+show_str if show_time else msg+show_str
    print('\r%s %d%%' %(show_str,percent*100),end='')

def setup_simulations(conditions_array, model):
    sim_array = []
    i = 0
    while (i < len(conditions_array)): #For all conditions
        properties = {} #create properties dictionary
        properties['temperature'] = float(conditions_array[i].temperature)
        properties['pressure'] = float(conditions_array[i].pressure)
        properties['equivalence_ratio'] = float(conditions_array[i].equi)
        properties['fuel'] = conditions_array[i].fuel
        properties['oxidizer'] = conditions_array[i].oxid

        sim_array.append(Simulation(i,properties,model)) #create simulation object and add it to the array
        i = i + 1

    return sim_array


############
# This function takes a set up simulation array, simulates it, and processes the results.  It returns a numpy array of the ignition delays
#
# sim_array: An array of set up simulation objects.  Comes in unsimulated and comes out simulated with results processed.
#
# ignition_delay: np array of the ignition delays from the setup_simulations.
#############


def simulate(sim_array):
        tau = [] #Ignition delays
        sample_points = [] #Sample information
        num = len(sim_array)
        for i,case in enumerate(sim_array): #Run simulations for original model and process results
            progress(float(i)/num,"Sim",1)
            tau.append(case.run_case())
            sample_points.append(case.process_results())

        ignition_delay = np.array(tau) #Turn tau array into a numpy array
        return ignition_delay
