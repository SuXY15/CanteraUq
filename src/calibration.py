from utils import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

pdiff_arr = 10**np.linspace(-4,0,10)

props['phi'] = 1.0
props['P'] = 10.0
props['T'] = 1000.

# = = = = = = = = = =
# calculator
def calculator(props_arr):
    NUM = len(props_arr)
    for i,props in enumerate(props_arr):
        print("Calculating data of c%d: %d/%d"%(rank,i,NUM))
        gas = load_mech(mech)
        pdiff = props['pdiff']
        name,mech = props['name'],props['mech']
        filename = comp_dir+name+"_calib.json"
        
        # brute force method for idt sensitivity
        sens_bf = np.zeros(gas.n_reactions)
        idt = get_ign(gas, props)
        for p in range(gas.n_reactions):
            gas.set_multiplier(1.+pdiff, p)
            pidt = get_ign(gas, props)
            sens_bf[p] = (np.log(pidt) - np.log(idt))/np.log(1.+pdiff)
            gas.set_multiplier(1.0, p)

        if gas.n_reactions<300:
            # approximate method for idt sensitivity
            idt,sens = get_ign_sens(gas, props, name='ap')
            p = np.argmax(abs(sens))
            gas.set_multiplier(1.+pdiff, p)
            pidt = get_ign(gas, props)
            gas.set_multiplier(1.0, p)
            tsens = (np.log(pidt) - np.log(idt))/np.log(1.+pdiff)
            sens_ap = sens*tsens/sens[p]
        else:
            sens_ap = sens_bf*0.99
        
        # saving
        props['sens_bf']= sens_bf.tolist()
        props['sens_ap']= sens_ap.tolist()
        json_writer(filename,props)

if __name__=="__main__":
    checkpath(comp_dir)
    checkpath(figs_dir)
    # = = = = = = = = = =
    # calculating
    props_arr = []
    for m,name,mech in mechs:
        filename = comp_dir+name+"_calib.json"
        if not checkexists(filename,delete=False):
            print("Data of %s not found, need %d calculation."%(name,len(pdiff_arr)))
            for pdiff in pdiff_arr:
                props['name'],props['mech'],props['pdiff'] = name, mech, pdiff
                props_arr.append(deepcopy(props))
        else:
            old_count = len(props_arr)
            props_old_arr = json_reader(filename)
            pdiff_old_arr = [props_old['pdiff'] for props_old in props_old_arr]
            for pdiff in pdiff_arr:
                if pdiff not in pdiff_old_arr:
                    props['name'],props['mech'],props['pdiff'] = name, mech, pdiff
                    props_arr.append(deepcopy(props))
            count = len(props_arr)-old_count
            print("Data of %s exists, need %d calculation."%(name,count))

    if len(props_arr):
        l,i = [int(li) for li in np.linspace(0, len(props_arr), num=size+1)],rank
        calculator(props_arr[l[i]:l[i+1]])
    
    if rank!=0: exit(0)

    # = = = = = = = = = =
    # showing results
    gas0 = load_mech(mechs[0][2])
    eqs0 = [r.equation for r in gas0.reactions()]
    fig = plt.figure("Difference Calibration",figsize=c2i(12,9))
    for m,name,mech in mechs:
        cprint("Reading data of %s"%mech,'b')
        gas = load_mech(mech)
        lr_arr,sens_list = [], []
        props_arr = json_reader(comp_dir+name+"_calib.json")
        for props in props_arr:
            pdiff = props['pdiff']
            sens_bf,sens_ap = np.array(props['sens_bf']),np.array(props['sens_ap'])

            # calculating length ratio
            norm_ap,norm_bf = normalize(sens_ap),normalize(sens_bf)
            lr_arr.append(np.dot(norm_ap,norm_bf))
            sens_list.append([norm_bf,norm_ap])
        
        sens_bf,sens_ap = sens_list[0]
        k = 10
        x = np.array(range(0,k))

        rank = sorted(range(len(sens_bf)), key=lambda k: -abs(sens_bf[k]))
        eqs = [r.equation for r in gas.reactions()]
        xtick = ['R%d'%pr for pr in parentRank(rank[:k], eqs, eqs0)]
        hbf = sens_bf[rank[:k]]
        hap = sens_ap[rank[:k]]
        
        plt.figure("Difference Calibration",figsize=c2i(12,9))
        # plt.title(r"Sensitivity of brute force method, $\delta$ calibration")
        plt.plot(pdiff_arr,1-np.array(lr_arr), color_arr[m]+'-')
        plt.xscale(r'log')
        plt.yscale(r'log')
        plt.xlabel(r'perturbation $\delta$')
        plt.ylabel(r'relative error')
        plt.legend(mech_arr,frameon=False)

        # plot
        fig2 = plt.figure("Sensitivity Calibration %s"%name,figsize=c2i(12,9))
        # plt.title("Sensitivity difference between methods, %s"%name)
        plt.bar(x-0.2,hbf,width=0.4)
        plt.bar(x+0.2,hap,width=0.4)
        plt.xticks(x,xtick,rotation=45)
        plt.xlabel(r'Reaction Index')
        plt.ylabel(r'Sensitivity')
        plt.legend([r'brute force',r'adjoint'],loc="upper right",frameon=False)
        save_figure(fig2, path=figs_dir+'calib_%s_sens.eps'%name)

    save_figure(fig, path=figs_dir+'calib_pdiff.eps')
    plt.show()

