import os, sys, time, json
import threading
import numpy as np
import cantera as ct
try:
    import matplotlib
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    font={'size':13}
    matplotlib.rc('font', **font)
except:
    "Note: matplotlib not found, are you on a server with no GUI ?"
from copy import deepcopy

# = = = = = = = = = =
# global settings
t0 = time.time()
threadLock = threading.Lock()
np.random.seed(0x532E532E53>>7)

def checkpath(path):
    '''
    Check if path or path of file exists, if not, create it
    '''
    if(path[-1]!='/'):
        path = os.path.split(path)[0]
    os.system("mkdir -p "+path)

def checkexists(file_name, delete=True):
    '''
    Check if file exists, if yes, delete it by default input
    '''
    if os.path.exists(file_name):
        if delete:
            try:
                os.remove(file_name)
            except:
                pass
        return True
    return False

def crosscheck(old_, new_):
    '''
    Check part of old_ not in new_ as useless
          part of new_ not in old_ as unfound
    '''
    useless, unfound = [], []
    for i in old_:
        if i not in new_:
            useless.append(i)
    for i in new_:
        if i not in old_:
            unfound.append(i)
    return useless, unfound

def cosine(a,b):
    return np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)

def parentRank(sRank,sList,pList):
    '''
    Find position of sList[sRank] in pList
    '''
    return [first(pList, lambda pi: pi==si)[0] for si in np.array(sList)[sRank]]

def childRank(pRank,sList,pList):
    '''
    Find position of pList[pRank] in sList
    '''
    return [first(sList, lambda si: si==pi)[0] for pi in np.array(pList)[pRank]]

def cprint(msg, color='r'):
    '''
    colorful print
    '''
    color_dict = {'r':31, 'g': 32, 'b': 34}
    print('\033[1;%dm%s\033[0m'%(color_dict[color],msg))

def progress(percent, msg='', width=40):
    '''
    Logging with progress bar
    '''
    percent = 1.0 if percent > 0.995 else percent + 0.005
    time_str = '|%6.1fm|'%((time.time()-t0)/60)
    show_str = ('[%%-%ds]' %width) %(int(width * percent)*"#")
    print('\r%s %d%%' %(msg+time_str + show_str, percent*100),end='')

def smooth(X,window_size=5):
    if len(X)<2*window_size:
        raise Exception("Length of data should be larger than 2 * window_size")
        sys.exit(0)
    N = len(X)
    return np.mean([X[i:N-window_size+i+1] for i in range(window_size)],axis=0)
    
def normalize(data):
    '''
    Normalize a vector
    '''
    return data/np.linalg.norm(data)

def cdiff(data):
    if len(data)==1: return np.array([0])
    return np.array([0]+list(np.diff(data)/2)) + np.array(list(np.diff(data)/2)+[0])

def diffMax(x,y):
    diff = cdiff(y)/cdiff(x)
    pos = np.argmax(diff)
    return pos, x[pos]

def curvMax(x,y,n=5,N=100):
    pos, x_val = diffMax(x, y)
    x, y = x[pos-n:pos+n], y[pos-n:pos+n]
    from scipy import interpolate
    x_new = np.linspace(x[0], x[-1], N)
    f = interpolate.interp1d(x, y, 'quadratic')
    return diffMax(x_new,f(x_new))

def first(the_iterable, condition = lambda x:True):
    '''
    Find first value satisfying the condition
    '''
    for i,li in enumerate(the_iterable):
        if condition(li):
            return i, li
    return None, None

def accum(the_iterable, condition = lambda x:x**2, sum_limit = 0.99):
    '''
    Find fist position satisfying the cumsum condition
    '''
    adder = 0.
    for i,li in enumerate(the_iterable):
        adder += condition(li)
        if adder > sum_limit:
            return i

def acquireLock(maxt, size, rank):
    '''
    Lock with time seperation between MPI threads
    '''
    while(int((time.time()-t0)/maxt)%size!=rank):
        time.sleep(1e-3)

def get_factors(uf):
    '''
    Get multiplier factor from uncertainty with lognormal
    '''
    return np.exp(np.random.normal(0.0, 1.0, len(uf))*np.log(uf)/3.)
    
def get_conditions(phi_arr, P_arr, T_arr):
    '''
    Get all conditions into one list
    '''
    condition_arr = []
    for phi in phi_arr:
        for P in P_arr:
            for T in T_arr:
                condition_arr.append([phi,P,T])
    return condition_arr

def get_props_arr(props, conditions):
    return [deepcopy(props) for (props['phi'],props['P'],props['T']) in conditions]

def get_mechs(mech_dir, mech_arr):
    '''
    Get all mechs into one list
    '''
    mech_list = []
    for i,mech_name in enumerate(mech_arr):
        mech_list.append([i,mech_name,mech_dir+mech_name+".cti"])
    return mech_list

def get_subspace(name, props, dim):
    subs_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d_as.json"%(name,
        UF,props['phi'],props['P'],props['T'],samplingSize)
    VT = np.array(json_reader(subs_name)[0]['subspace'])
    return np.transpose(VT[:dim])  # d x s

def get_Xy(name, props, uf):
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
        UF,props['phi'],props['P'],props['T'],samplingSize)
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]

    idt = np.array([props['idt'] for props in props_arr])
    idx = [i for i,iidt in enumerate(idt) if iidt>0]
    f = np.transpose([np.log10(idt[idx]),])
    X_raw = [3.*np.log(props_arr[i]['factor'])/np.log(uf) for i in idx]
    return X_raw, f # n x d, n x 1

def get_sub_plots(num):
    fig, axs = plt.subplots(len(phi_arr), len(P_arr),figsize=(12,9),num=num)
    fig.subplots_adjust(left=0.05,bottom=0.05,top=0.95,right=0.95,hspace=0.,wspace=0.)
    fig.canvas.set_window_title(num)
    fig.suptitle(num)
    return fig,axs

def set_sub_plots(axs,xlabel,ylabel,legend,xlim=None,ylim=None,xscale=None,yscale=None,):
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            if(p!=0):               axs[i,p].get_yaxis().set_visible(False)
            if(i!=len(phi_arr)-1):  axs[i,p].get_xaxis().set_visible(False)
            if i==2 and p==1:       axs[i,p].set_xlabel(xlabel)
            if i==1 and p==0:       axs[i,p].set_ylabel(ylabel)
            if i==1 and p==1:       axs[i,p].legend(legend,loc="lower right",frameon=False)
            if xlim:    axs[i,p].set_xlim(xlim)
            if ylim:    axs[i,p].set_ylim(ylim)
            if xscale:  axs[i,p].set_xscale(xscale)
            if yscale:  axs[i,p].set_yscale(yscale)
            axs[i,p].set_title(r"$\phi=%.1f, P=%.0fatm$"%(phi,P), pad=-20)

def save_figure(fig, path=None):
    fig.savefig(path, format='png', bbox_inches='tight', transparent=True, dpi=600)

def sub_mech_writer(original_mech, species, reactions, name=''):
    '''
    Save sub mechanism to file with bool index
    Eg: sub_mech_writer(mech, species, reactions, '_CH4')
    '''
    # prepare
    sub_mech_name = original_mech + '.' + name
    all_species = ct.Species.listFromFile(original_mech)
    all_reactions = ct.Reaction.listFromFile(original_mech)

    # generate index
    species_name = [s.name for s in species]
    reactions_equation = [r.equation for r in reactions]
    species_index = [1 if s.name in species_name else 0 for s in all_species]
    reactions_index = [1 if r.equation in reactions_equation else 0 for r in all_reactions]

    # save index
    sub_mech = open(sub_mech_name, 'w')
    sub_mech.writelines(["%s "%si  for si in species_index])
    sub_mech.writelines("\n")
    sub_mech.writelines(["%s "%ri  for ri in reactions_index])

def sub_mech_reader(original_mech, name=''):
    '''
    Read sub mechanism from file with bool index
    Eg: species,reactions = sub_mech_reader(mech, '_CH4')
    '''
    # prepare
    sub_mech_name = original_mech + '.' + name
    all_species = ct.Species.listFromFile(original_mech)
    all_reactions = ct.Reaction.listFromFile(original_mech)
    
    # get index
    sub_mech = open(sub_mech_name)
    species_index = [int(i) for i in sub_mech.readline().split()]
    reactions_index = [int(i) for i in sub_mech.readline().split()]

    # get species and reactions
    species = [s for i,s in enumerate(all_species) if species_index[i]]
    reactions = [r for i,r in enumerate(all_reactions) if reactions_index[i]]

    return species, reactions

def load_mech(mech_path='gri30.cti'):
    '''
    Recognize mech and sub mech by its name and load mechanism
    '''
    mech_split = mech_path.split('.')
    suffix = mech_split[-1]
    if suffix in ['xml', 'cti']:
        species = ct.Species.listFromFile(mech_path)
        reactions = ct.Reaction.listFromFile(mech_path)
    else:
        original_mech = mech_split[0]+'.'+mech_split[1]
        species, reactions = sub_mech_reader(original_mech, mech_split[-1])
    return ct.Solution(thermo='IdealGas', kinetics='GasKinetics', species=species, reactions=reactions)

def load_uncertainty(fileName, UF=2):
    '''
    Load uncertainty factor to a list
    '''
    uncertainty_coefs = []
    if fileName.split('.')[-1]=="txt":
        with open(fileName,'r') as f:
            for line in f.readlines():
                uncertainty_coefs.append(float(line.split()[1]))
    else:
        gas = load_mech(fileName)
        uncertainty_coefs = [UF for i in range(gas.n_reactions)]
    return uncertainty_coefs

def json_writer(filename, data, fmt=False):
    '''
    Write sensitivity data as json format to one line of file
    '''
    threadLock.acquire()
    with open(filename, "a+") as f:
        f.writelines(json.dumps(data)+"\r\n")
    threadLock.release()

def json_reader(filename, fmt=False):
    '''
    Read sensitivity data file, every line convert to json format
    '''
    data = []
    with open(filename, "r") as f:
        for line in f.readlines():
            data.append(json.loads(line))
        return data

def config_writer(config_file, data):
    '''
    Save config
    '''
    threadLock.acquire()
    with open(config_file, "w") as f:
        f.writelines(json.dumps(data, indent=4, separators=(',', ':'))+"\r\n")
    threadLock.release()

def config_reader(config_file):
    '''
    Read config
    '''
    with open(config_file, "r") as f:
        json_str = ""
        for line in f.readlines():
            json_str += line
    return json.loads(json_str)


def conditions_writer(mech_dir,props,conditions_list,note=''):
    filename = mech_dir+note+"conditions.txt"
    with open(filename, "w") as f:
        for props['phi'],props['P'],props['T'] in conditions_list:
            f.writelines("CONV"+"\r\n")
            f.writelines("PRES %.1f"%props['P']+"\r\n")
            f.writelines("TEMP %.1f"%props['T']+"\r\n")
            f.writelines("TIME %.1f"%props['tot']+"\r\n")
            f.writelines("DELT 1.E-4"+"\r\n")
            for fuel in props['fuel'].split(','):
                f.writelines("FUEL "+fuel.replace(':',' ')+"\r\n")
            for oxid in props['air'].split(','):
                f.writelines("OXID "+oxid.replace(':',' ')+"\r\n")
            f.writelines("EQUI %.1f"%props['phi']+"\r\n")
            f.writelines("END\r\n"+"\r\n")

def gas_reactor(gas, props, name=''):
    gas.set_equivalence_ratio(props['phi'], props['fuel'], props['air'])
    gas.TP = props['T'], props['P']*ct.one_atm
    return ct.IdealGasReactor(gas, name=name+'IGR') if props['type'] == 'UV' \
        else ct.IdealGasConstPressureReactor(gas, name=name+'IGCPR')

def get_ign(gas, props, name=''):
    '''
    Get iginition delay time
    '''
    r = gas_reactor(gas, props, name)

    # simulation
    sim = ct.ReactorNet([r])

    # set the tolerances for the solution and for the sensitivity coefficients
    if props['T'] >= 1200:
        sim.rtol = 1.0e-12
        sim.atol = 1.0e-15

    states = ct.SolutionArray(gas, extra=['t'])
    try:
        while sim.time < props['tot']:
            sim.step()
            states.append( r.thermo.state, t=sim.time)
    except:
        cprint("Error Occured In Cantera Solver!")
        return 0
    
    # handle results
    if max(states.T) < props['T'] + 200:
        cprint("Not ignited. Highest T:{:.2f}".format(max(states.T)))
        return 0
    else:
        ip, ign_delay = curvMax(states.t, states.T)
        return ign_delay

def getTcurv(gas, props, name=''):
    '''
    Get Temperature curve during an auto-ignition
    '''
    r = gas_reactor(gas, props, name)

    # simulation
    sim = ct.ReactorNet([r])
    if props['T'] >= 1200:
        sim.rtol = 1.0e-9
        sim.atol = 1.0e-12

    states = ct.SolutionArray(gas, extra=['t'])
    try:
        while sim.time < props['tot']:
            sim.step()
            states.append( r.thermo.state, t=sim.time)
    except:
        cprint("Error Occured In Cantera Solver!")
        return 0
    data = {'t':states.t, 'T':states.T}
    return data

def get_ign_sens(gas, props, name=''):
    '''
    Get iginition delay time and sensitivity
    '''
    r = gas_reactor(gas, props, name)

    sim = ct.ReactorNet([r])
    for i in range(gas.n_reactions):
        r.add_sensitivity_reaction(i)

    # set the tolerances for the solution and for the sensitivity coefficients
    # if props['T']>=1200:
    #     sim.rtol = 1.0e-9
    #     sim.atol = 1.0e-12
    #     sim.rtol_sensitivity = 1.0e-9
    #     sim.atol_sensitivity = 1.0e-12

    sens = []
    states = ct.SolutionArray(gas, extra=['t'])
    try:
        sim.verbose = False
        while sim.time < props['tot']:
            sim.step()
            sens.append(sim.sensitivities()[2,:])
            states.append(r.thermo.state, t=sim.time)
    except:
        cprint("Error Occured In Cantera Solver!")
        return 0, np.zeros(gas.n_reactions)
    sens = np.array(sens)
    
    # handle results
    if max(states.T) < props['T'] + 200:
        cprint("Not ignited.\nHighest Temperature:{:.2f}".format(max(states.T)))
        return 0, np.zeros(gas.n_reactions)
    else:
        dT = np.diff(states.T)/np.diff(states.t)
        sens = np.array(sens)
        ign_pos = np.argmax(dT)
        ign_delay = states.t[ign_pos]
        return ign_delay, sens[ign_pos,:]

def get_ign_sens_bf(gas, props, name=''):
    '''
    Get iginition delay time and sensitivity by brute force method
    '''
    sens = list(np.zeros(gas.n_reactions))
    idt = get_ign(gas,props)
    for p in range(gas.n_reactions):
        gas.set_multiplier(props['factor'][p]*(1+props['pdiff']), p)
        pidt = get_ign(gas, props)
        sens[p] = (np.log(pidt) - np.log(idt))/np.log(1+props['pdiff'])
        gas.set_multiplier(props['factor'][p], p)
    return idt,sens

def get_fs(gas, props, name=''):
    '''
    Get free flame speed
    '''
    gas.set_equivalence_ratio( props['phi'], props['fuel'], props['air'] )
    gas.TP = props['T'], props['P']*ct.one_atm
    
    sim = ct.FreeFlame(gas, width=props['width'])
    sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)

    # simulation
    sim.transport_model = 'Multi'
    sim.solve(loglevel=0, auto=True)
    
    flame_speed = sim.u[0]
    return flame_speed

def get_fs_sens(gas, props, name=''):
    '''
    Get free flame speed and sensitivity
    '''
    gas.set_equivalence_ratio( props['phi'], props['fuel'], props['air'] )
    gas.TP = props['T'], props['P']*ct.one_atm
    
    # reaction setting
    sim = ct.FreeFlame(gas, width=props['width'])
    sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)

    # simulation
    sim.transport_model = 'Multi'
    sim.solve(loglevel=0, auto=True)
    
    # handle results
    sens = sim.get_flame_speed_reaction_sensitivities()
    flame_speed = sim.u[0]
    return flame_speed, sens

# plt.ion()
# fig = plt.figure()
# ax = plt.subplot()

    # ax.plot(sens,'b')
    # plt.draw()
    # plt.pause(1e-17)
    # ax.lines.pop(0)

color_arr = ['k','b','r','c']
symbol_arr = ['s','o','d','v','^']
config_file = "nothing.config"
try:
    config_file = sys.argv[1]+".config"
    config_dict = config_reader(config_file)
    # loading configs
    UF = config_dict['UF']
    maxt = config_dict['maxt']
    pdiff = config_dict['pdiff']
    samplingSize = config_dict['samplingSize']
    comp_dir = config_dict['comp_dir']
    sens_dir = config_dict['sens_dir']
    acts_dir = config_dict['acts_dir']
    resp_dir = config_dict['resp_dir']
    figs_dir = config_dict['figs_dir']
    mech_dir = config_dict['mech_dir']
    mech_arr = config_dict['mech_arr']
    phi_arr  = config_dict['phi_arr']
    P_arr = config_dict['P_arr']
    T_arr = config_dict['T_arr']
    props = config_dict['props']
    mechs = get_mechs(mech_dir,mech_arr)
    conditions = get_conditions(phi_arr,P_arr,T_arr)
    cprint("%s Loaded"%config_file,'g')
except:
    cprint("Loading %s failed."%config_file)
