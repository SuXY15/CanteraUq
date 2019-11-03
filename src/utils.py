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
    matplotlib.rc('text', usetex=True)
    # matplotlib.rc('ps', usedistiller='ghostscript')
except:
    "Note: matplotlib not found, are you on a server with no GUI?"
from copy import deepcopy

# = = = = = = = = = =
# global settings
t0 = time.time()
threadLock = threading.Lock()
np.random.seed(0x532E532E53>>7)
line_arr = ('-','--','-.',':')
color_arr = ('k','r','b','m')
symbol_arr = ('s','o','v','^')

# = = = = = = = = = =
# message in cli
def cprint(msg, color='r'):
    """ Colorful print
    """
    color_dict = {'r':31, 'g': 32, 'b': 34}
    print('\033[1;%dm%s\033[0m'%(color_dict[color],msg))

def progress(percent, msg='', width=40):
    """ Logging with progress bar
    """
    percent = 1.0 if percent > 0.995 else percent + 0.005
    time_str = '|%6.1fm|'%((time.time()-t0)/60)
    show_str = ('[%%-%ds]' %width) %(int(width * percent)*"#")
    print('\r%s %d%%' %(msg+time_str + show_str, percent*100),end='')


# = = = = = = = = = =
# for directories and files
def acquireLock(maxt, size, rank):
    """ Lock with time seperation between MPI threads
    """
    while(int((time.time()-t0)/maxt)%size!=rank):
        time.sleep(1e-3)

def checkpath(path):
    """ Check if path or path of file exists, if not, create it
    """
    path = os.path.split(path)[0] if path[-1]!='/' else path
    os.system("mkdir -p "+path)

def checkexists(filename, delete=True):
    """ Check if file exists, if yes, delete it by default input
    """
    if os.path.exists(filename):
        if delete:
            try:
                os.remove(filename)
            except:
                pass
        return True
    return False

def json_writer(filename, data, fmt=False):
    """ Write data as json format to one line of file
    """
    threadLock.acquire()
    with open(filename, "a+") as f:
        f.writelines(json.dumps(data)+"\r\n")
    threadLock.release()

def json_reader(filename, fmt=False):
    """ Read json data file, every line convert to json format
    """
    data = []
    with open(filename, "r") as f:
        for line in f.readlines():
            data.append(json.loads(line))
        return data

def config_writer(config_file, data):
    """ Save config
    """
    threadLock.acquire()
    with open(config_file, "w") as f:
        f.writelines(json.dumps(data, indent=4, separators=(',', ':'))+"\r\n")
    threadLock.release()

def config_reader(config_file):
    """ Read config
    """
    json_str = ""
    with open(config_file, "r") as f:
        for line in f.readlines():
            json_str += line
    return json.loads(json_str)

# = = = = = = = = = =
# find conditional rank
def first(iterable, condition = lambda x: True):
    """ Find first value satisfying the condition
    """
    for i,li in enumerate(iterable):
        if condition(li):
            return i, li
    return None, None

def accum(iterable, condition = lambda x: x**2, sum_limit = 0.99):
    """ Find fist position satisfying the cumsum condition
    """
    adder = 0.
    for i,li in enumerate(iterable):
        adder += condition(li)
        if adder > sum_limit:
            return i

def parentRank(sRank, sList, pList):
    """ Find position of sList[sRank] in pList
    """
    return [first(pList, lambda pi: pi==si)[0] for si in np.array(sList)[sRank]]

def childRank(pRank, sList, pList):
    """ Find position of pList[pRank] in sList
    """
    return [first(sList, lambda si: si==pi)[0] for pi in np.array(pList)[pRank]]

def crosscheck(old_arr, new_arr):
    """ Check part of old_arr not in new_arr as useless
              part of new_arr not in old_arr as unfound
    """
    useless, unfound = [], []
    for i in old_arr:
        if i not in new_arr:
            useless.append(i)
    for i in new_arr:
        if i not in old_arr:
            unfound.append(i)
    return useless, unfound

# = = = = = = = = = =
# algebra
def normalize(data):
    """ Normalize a vector
    """
    return data/np.linalg.norm(data)

def cosine(a, b):
    """ Calculate consine of two vector
    """
    return np.dot(normalize(a), normalize(b))

def cdiff(data):
    """ Central differential of data
    """
    if len(data)==1: return np.array([0])
    d = np.diff(data)
    return np.array([d[0]] + list((d[1:]+d[:-1])/2) + [d[-1]])

def diffMax(x, y):
    """ Find max position of dy/dx
    """
    pos = np.argmax(cdiff(y)/cdiff(x))
    return pos, x[pos]

def curvMax(x, y, n=5, N=100):
    """ Find max position of dy/dx with high resolution
    """
    from scipy import interpolate
    l = len(x)-1
    pos, x_val = diffMax(x, y)
    lpos = 0 if pos-n<0 else pos-n
    rpos = l if pos+n>l else pos+n
    x, y = x[lpos:rpos], y[lpos:rpos]
    x_new = np.linspace(x[0], x[-1], N)
    f = interpolate.interp1d(x, y, 'quadratic')
    return diffMax(x_new,f(x_new))

# = = = = = = = = = =
# for graphical / saving figures
def c2i(*tupl):
    """ Convert centermeter to inch
    """
    return tuple(i/2.54 for i in tupl)

def get_sub_plots(num):
    """ Get len(phi_arr)xlen(P_arr) subplots, 3x3 recommended
    """
    fig, axs = plt.subplots(len(P_arr), len(phi_arr), figsize=c2i(25.8,15), num=num)
    fig.subplots_adjust(left=0.05,bottom=0.05,top=0.95,right=0.95,hspace=0.,wspace=0.)
    fig.canvas.set_window_title(num)
    return fig, axs

def set_sub_plots(axs, xlabel, ylabel, legend, xlim=None, ylim=None, xscale=None, yscale=None,):
    """ Set subplots
    """
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            if(p!=0):              axs[i,p].get_yaxis().set_visible(False)
            if(i!=len(phi_arr)-1): axs[i,p].get_xaxis().set_visible(False)
            if i==2 and p==1:      axs[i,p].set_xlabel(xlabel)
            if i==1 and p==0:      axs[i,p].set_ylabel(ylabel)
            if i==1 and p==1:      axs[i,p].legend(legend, loc="lower right", frameon=False)
            if xlim:   axs[i,p].set_xlim(xlim)
            if ylim:   axs[i,p].set_ylim(ylim)
            if xscale: axs[i,p].set_xscale(xscale)
            if yscale: axs[i,p].set_yscale(yscale)
            axs[p,i].set_title(r"$\phi=%.1f, P=%.0fatm$"%(phi,P), pad=-20)

def save_figure(fig, path=None):
    """ Save figure with handle
    """
    fig.savefig(path, bbox_inches='tight', transparent=True, dpi=300)

# = = = = = = = = = =
# numerous functions for convenience
def get_factors(UF):
    """ Get multiplier factor from uncertainty with lognormal
    """
    return np.exp(np.random.normal(0.0, 1.0, len(UF))*np.log(UF)/3.)
    
def get_conditions(phi_arr, P_arr, T_arr):
    """ Get all conditions into one list
    """
    condition_arr = []
    for phi in phi_arr:
        for P in P_arr:
            for T in T_arr:
                condition_arr.append([phi,P,T])
    return condition_arr

def get_props_arr(props, conditions):
    """ Get array of props
    """
    return [deepcopy(props) for (props['phi'],props['P'],props['T']) in conditions]

def get_mechs(mech_dir, mech_arr):
    """ Get all mechs into one list
    """
    mech_list = []
    for i,mech_name in enumerate(mech_arr):
        mech_list.append([i,mech_name,mech_dir+mech_name+".cti"])
    return mech_list

def get_subspace(name, props, dim=1):
    """ Get subspace data, dim is required
    """
    filename = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d_as.json"%(name,
                    props['UF'],props['phi'],props['P'],props['T'],props['S'])
    VT = np.array(json_reader(filename)[0]['subspace'])
    return np.transpose(VT[:dim])  # shape(d, s)

def get_Xy(name, props, uf):
    """ Get input X and output y
    """
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                    props['UF'],props['phi'],props['P'],props['T'],props['S'])
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]
    idt = np.array([props['idt'] for props in props_arr])
    idx = [i for i,iidt in enumerate(idt) if iidt>0]
    f = np.transpose([np.log10(idt[idx[:360]]),])
    X_raw = [3.*np.log(props_arr[i]['factor'])/np.log(uf) for i in idx[:360]]
    return X_raw, f # shape(n, d), shape(n, 1)

def load_mech(mech='gri30.cti'):
    """ Load mechanism
    """
    return ct.Solution(mech)

def load_uncertainty(filename, UF=2):
    """ Load uncertainty factor to a list, if not provided, input a mech
    """
    UFs = []
    if filename.split('.')[-1]=="txt":
        with open(filename,'r') as f:
            for line in f.readlines():
                UFs.append(float(line.split()[1]))
    else:
        gas = load_mech(filename)
        UFs = [UF for i in range(gas.n_reactions)]
    return UFs

def conditions_writer(mech_dir,props,conditions_list):
    """ write conditions file for DRG
    """
    filename = mech_dir+"conditions.txt"
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

# = = = = = = = = = =
# Cantera simulator
def gas_reactor(gas, props, name=''):
    """ Get a GasReactor
    """
    gas.set_equivalence_ratio(props['phi'], props['fuel'], props['air'])
    gas.TP = props['T'], props['P']*ct.one_atm
    return ct.IdealGasReactor(gas, name=name+'IGR') if props['type'] == 'UV' \
            else ct.IdealGasConstPressureReactor(gas, name=name+'IGCPR')

def getTcurv(gas, props, name=''):
    """
    Get Temperature curve during an auto-ignition
    """
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

def get_ign(gas, props, name=''):
    """ Get iginition delay time
    """
    data = getTcurv(gas, props, name)
    maxT =  max(data['T'])
    if maxT < props['T'] + 200:
        cprint("Not ignited. Highest T:{:.2f}".format(maxT))
        return 0
    else:
        _, idt = curvMax(data['t'], data['T'])
        return idt

def get_ign_sens(gas, props, name=''):
    """ Get iginition delay time and sensitivity
    """
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
        cprint("Not ignited. Highest T:{:.2f}".format(max(states.T)))
        return 0, np.zeros(gas.n_reactions)
    else:
        pos, idt = diffMax(states.t, states.T)
        return idt, sens[pos,:]

def get_ign_sens_bf(gas, props):
    """ Get iginition delay time and sensitivity by brute force method
    """
    sens = list(np.zeros(gas.n_reactions))
    idt = get_ign(gas, props)
    for p in range(gas.n_reactions):
        gas.set_multiplier(props['factor'][p]*(1+props['pdiff']), p)
        pidt = get_ign(gas, props)
        sens[p] = (np.log(pidt) - np.log(idt))/np.log(1+props['pdiff'])
        gas.set_multiplier(props['factor'][p], p)
    return idt,sens

def get_fs(gas, props, name=''):
    """ Get free flame speed
    """
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
    """
    Get free flame speed and sensitivity
    """
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

def set_factors(gas, factors):
    """
    Set factors for cantera gas simulators
    """
    gas.set_multiplier(1.0)
    for i,factor in enumerate(factors):
        gas.set_multiplier(factor, i)
    return gas

# loading configs
if len(sys.argv)>1:
    config_file = 'data/'+sys.argv[1]+".config"
    try:
        config_dict = config_reader(config_file)
        locals().update(config_dict)
        mechs = get_mechs(mech_dir,mech_arr)
        conditions = get_conditions(phi_arr,P_arr,T_arr)
        cprint("%s loaded"%config_file,'g')
        props['UF'] = UF
        props['S'] = samplingSize
    except:
        cprint("Load %s failed, please check."%config_file,'r')
