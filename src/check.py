from utils import *

check = [int(i) for i in "0 0 1 1 1 1 1 1".split()]

# check all paths
def check_path():
    print("\nCheck all paths:")
    print("\t%s"%mech_dir); checkpath(mech_dir)
    print("\t%s"%comp_dir); checkpath(comp_dir)
    print("\t%s"%sens_dir); checkpath(sens_dir)
    print("\t%s"%acts_dir); checkpath(acts_dir)
    print("\t%s"%resp_dir); checkpath(resp_dir)

# check meshs
def check_mech():
    print("\nCheck if all mechs ready:")
    for i,name,mech in mechs:
        print("\t%s in %s"%(name,mech),end='')
        try:
            load_mech(mech)
            print(" ready.")
        except:
            cprint(" error!")

# check compare
def check_comp(file_name, conditions, delete=False):
    data_arr = json_reader(file_name)
    props_arr = [d['props'] for d in data_arr]
    old_conditions = [[props['phi'],props['P'],props['T']] for props in props_arr]
    useless,unfound = crosscheck(old_conditions, conditions)
    for c in useless:
        print("\tUseless: phi %s P %s T %s"%(c[0],c[1],c[2]))
    new_props_arr = []
    new_props = props_arr[0]
    for c in unfound:
        new_props['phi'],new_props['P'],new_props['T'] = c
        new_props_arr.append(deepcopy(new_props))
    if delete:
        checkexists(file_name,delete=True)
        for props in props_arr:
            if [props['phi'],props['P'],props['T']] in conditions:
                json_writer(file_name,{'props':props,'tdata':[]})
            
    return new_props_arr

# check sensitivity
def check_sens(name,mech,delete=False):
    file_name = sens_dir+name+".json"
    data_arr = json_reader(file_name)
    props_arr = [d['props'] for d in data_arr]
    old_conditions = [[props['phi'],props['P'],props['T']] for props in props_arr]
    useless,unfound = crosscheck(old_conditions, conditions)
    for c in useless:
        print("\tUseless: phi %s P %s T %s"%(c[0],c[1],c[2]))
    new_props_arr = []
    new_props = props_arr[0]
    for c in unfound:
        new_props['phi'],new_props['P'],new_props['T'] = c
        new_props_arr.append(deepcopy(new_props))
    if delete:
        checkexists(file_name,delete=True)
        for d in data_arr:
            props, tdata = d['props'], d['tdata']
            if [props['phi'],props['P'],props['T']] in conditions:
                json_writer(file_name,{'props':props,'tdata':tdata})
    return new_props_arr

# check active subspace
def check_acts(delete=False):
    old_files = []
    for fi in os.listdir(acts_dir):
        if len(fi)>5 and fi[-7:-5] != 'as':
            old_files.append(fi)
    sup_files = []
    for i,name,mech in mechs:
        for c in conditions:
            file_name = "%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                        UF,c[0],c[1],c[2],samplingSize)
            sup_files.append(file_name)
    useless, unfound = crosscheck(old_files, sup_files)
    for fi in useless:
        print("  Useless: %s"%fi)
        checkexists(file_name,delete=delete)
    for fi in unfound:
        print("  Unfound: %s"%fi)
    return unfound

def check_continuity():
    m, name, mech = mechs[3]
    gas = load_mech(mech)
    props['factor'] = list(np.ones(gas.n_reactions))
    props['pdiff'] = pdiff
    props['phi'], props['P'], props['T'] = 0.5, 1.0, 1600
    
    T1 = getTcurv(gas, props)
    gas.set_multiplier(1.+props['pdiff'], 0)
    T2 = getTcurv(gas, props)

    ip1,idt1 = diffMax(T1['t'], T1['T'])
    ip2,idt2 = diffMax(T2['t'], T2['T'])
    sens = (np.log(idt2) - np.log(idt1))/np.log(1+props['pdiff'])
    print(sens)

    id1 = range(ip1-5, ip1+5)
    id2 = range(ip2-5, ip2+5)
    fig = plt.figure("High Temperature IDT",figsize=(6,4.5))
    plt.title("High Temperature IDT")
    plt.plot(T1['t'][id1]*1e6,(cdiff(T1['T'])/cdiff(T1['t']))[id1],'b')
    plt.plot(T2['t'][id2]*1e6,(cdiff(T2['T'])/cdiff(T2['t']))[id2],'r')

    ip1,idt1 = curvMax(T1['t'], T1['T'])
    ip2,idt2 = curvMax(T2['t'], T2['T'])
    sens = (np.log(idt2) - np.log(idt1))/np.log(1+props['pdiff'])
    print(sens)

    plt.xlabel(r'$t / \mu s$')
    plt.ylabel(r'$\frac{\partial T}{\partial t}$')
    plt.legend([r'original',r'perturbated'], frameon=False)
    save_figure(fig, path=figs_dir+"sens_highTcomp.png")
    plt.show()

def checkTcurv():
    m, name, mech = mechs[1]
    props = config_dict['props']

    props['phi'], props['P'], props['T'] = 1.0, 1.0, 800.0
    lineIdx = 69

    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]
    
    plt.figure()
    for tdata in tdata_arr:
        plt.plot(tdata)

    tdata, props = tdata_arr[lineIdx], props_arr[lineIdx]
    props['tot']=10.0

    gas.set_multiplier(1.0)
    for j,factor in enumerate(props['factor']):
        gas.set_multiplier(factor, j)
    sens = np.zeros(gas.n_reactions)
    idt = get_ign(gas, props)

    plt.figure()
    for p in props['mri']:
        data = getTcurv(gas, props)
        plt.plot(data['t'], data['T'])
        plt.plot(data['t'], cdiff(data['T'])/cdiff(data['t']))
    
    # plt.figure()
    # for p in props['mri']:
    #     gas.set_multiplier(props['factor'][p]*(1+pdiff), p)
    #     pidt = get_ign(gas, props)
    #     sens[p] = (np.log(pidt) - np.log(idt))/np.log(1+pdiff)
    #     print("%6.3e %6.3e"%(idt, pidt))
    #     gas.set_multiplier(props['factor'][p], p)

    # plt.plot(sens)
    plt.show()

def revise_sampling():
    m, name, mech = mechs[1]
    props = config_dict['props']
    props['phi'], props['P'], props['T'] = 1.0, 1.0, 700.0
    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    sens_data_arr = json_reader(file_name)
    
    new_data_arr = []
    checkexists(file_name, delete=True)
    for sd in sens_data_arr:
        tdata, props = sd['tdata'], sd['props']
        
        print(np.sum(np.abs(tdata)))

        if np.sum(np.abs(tdata))>100:
            props['tot']=10.0
            gas.set_multiplier(1.0)
            for j,factor in enumerate(props['factor']):
                gas.set_multiplier(factor, j)
            sens = np.zeros(gas.n_reactions)
            
            idt = get_ign(gas, props)
            for p in props['mri']:
                gas.set_multiplier(props['factor'][p]*(1+pdiff), p)
                pidt = get_ign(gas, props)
                sens[p] = (np.log(pidt) - np.log(idt))/np.log(1+pdiff)
                gas.set_multiplier(props['factor'][p], p)

            json_writer(file_name, {'props': deepcopy(props),'tdata': sens.tolist()})
        else:
            json_writer(file_name, sd)


if __name__=="__main__":
    # # check compare
    # for i,name,mech in mechs:
    #     print("\nCheck if %s compare data:"%name)
    #     for props in check_comp(name,mech,True):
    #         print("\tUnfound: phi %s P %s T %s"%(props['phi'],props['P'],props['T']))

    # # check sensitivity
    # for i,name,mech in mechs:
    #     print("\nCheck if %s sensitivity data:"%name)
    #     for props in check_sens(name,mech,True):
    #         print("\tUnfound: phi %s P %s T %s"%(props['phi'],props['P'],props['T']))

    # # check active subspace
    # check_acts()

    # # check high temperature sensitivity
    # check_continuity()

    # check sensitivity cut off by 1.0s
    checkTcurv()

    # # revise the error from cut off at 1.0s
    # revise_sampling()