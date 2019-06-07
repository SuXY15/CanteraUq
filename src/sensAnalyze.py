from utils import *
from check import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# = = = = = = = = = =
# calculator
def calculator(id, calc_arr):
    sens_dict = {'props':[], 'tdata':[]}
    NUM = len(calc_arr)
    for i,[name,mech,props] in enumerate(calc_arr):
        file_name = sens_dir+name+".json"
        gas = load_mech(mech)
        print("|{0:07.3f}s| calculator {1}: {2}/{3}...".format(time.time()-t0, id, i, NUM))
        # brute force method for idt sensitivity
        props['idt'],tdata = get_ign_sens_bf(gas,props)
    # if gas.n_reactions<300:
        # approximate method for idt sensitivity
        idt,sens = get_ign_sens(gas,props,name='ap')
        p = np.argmax(abs(sens))
        gas.set_multiplier(1.+props['pdiff'], p)
        pidt = get_ign(gas, props)
        gas.set_multiplier(1.0, p)
        tsens = (np.log(pidt) - np.log(idt))/np.log(1.+props['pdiff'])
        sdata = sens*tsens/sens[p]
    # else:
    #     sdata = tdata*0.99
        # save data
        if props['idt'] != 0:
            sens_dict['props'] = deepcopy(props)
            sens_dict['tdata'] = deepcopy(tdata)
            sens_dict['sdata'] = deepcopy(sdata.tolist())
            acquireLock(maxt, size, rank)
            json_writer(file_name, sens_dict)

if __name__=="__main__":
    checkpath(sens_dir)
    checkpath(figs_dir+"sens_comp/")
    # = = = = = = = = = =
    # calculating
    calc_arr = []
    for m,name,mech in mechs:
        gas = load_mech(mech)
        file_name = sens_dir+name+".json"
        if not checkexists(file_name, delete=False):
            props['factor'] = list(np.ones(gas.n_reactions))
            props['pdiff'] = pdiff
            for props in get_props_arr(props, conditions):
                calc_arr.append([deepcopy(name),deepcopy(mech),deepcopy(props)])
        else:
            for props in check_sens(name,mech):
                calc_arr.append([deepcopy(name),deepcopy(mech),deepcopy(props)])

    if len(calc_arr):
        l,i = [int(li) for li in np.linspace(0, len(calc_arr), num=size+1)],rank
        calculator('c%d'%i, calc_arr[l[i]:l[i+1]])

    if rank!= 0: sys.exit(0)

    # = = = = = = = = = =
    # calculate global main reactions
    gas0 = load_mech(mechs[0][2])
    eqs0 = [r.equation for r in gas0.reactions()]
    main_reaction = set()
    for m,name,mech in mechs:
        gas = load_mech(mech)
        filename = sens_dir+name+".json"
        main_reaction_m = set()
        eqs = [r.equation for r in gas.reactions()]
        for data in json_reader(filename):
            props = data['props']
            tdata = normalize(data['tdata'])
            sdata = normalize(data['sdata'])
            # accum for main reactions
            rank = np.argsort(-abs(sdata))
            if np.sum(np.abs(sdata))>0.5:
                pos = accum(sdata[rank], sum_limit=0.99 if props['T']<1200 else 0.5)
                pRank = parentRank(rank[0:pos],eqs,eqs0)
                main_reaction_m.update(pRank)
            print("%s %.1f, %.0fatm, %.0fK: %d %d"%
                (name, props['phi'], props['P'], props['T'], len(main_reaction_m),len(main_reaction)))
        main_reaction.update(main_reaction_m)
    main_reaction = sorted(list(main_reaction))

    # = = = = = = = = = =
    # show results
    for m,name,mech in mechs:
        cprint("Loading %s"%name,'b')
        gas = load_mech(mech)
        eqs = [r.equation for r in gas.reactions()]
        rank = np.array(childRank(main_reaction, eqs, eqs0))
        idx = [i for i,r in enumerate(rank) if r is not None]
        rank = rank[idx].astype('int')
        x = np.array(range(0,len(main_reaction)))[idx]
        xtick = np.array(['R%d'%pr for pr in main_reaction])[idx]

        checkexists(sens_dir+name+"_mr.json",delete=True)
        json_writer(sens_dir+name+"_mr.json",{"mr":rank.tolist()})

        data_arr = json_reader(sens_dir+name+".json")
        props_arr = [d['props'] for d in data_arr]
        tdata_arr = np.array([normalize(d['tdata']) for d in data_arr])
        sdata_arr = np.array([normalize(d['sdata']) for d in data_arr])
        err_arr = []
        for t,T in enumerate(T_arr):
            fig, AX = get_sub_plots("Sensitivity Comparation, %s, T=%.0fK"%(name,T))
            props_tarr = [props_arr[i] for i,p in enumerate(props_arr) if p['T'] == T]
            tdata_tarr = [tdata_arr[i] for i,p in enumerate(props_arr) if p['T'] == T]
            sdata_tarr = [sdata_arr[i] for i,p in enumerate(props_arr) if p['T'] == T]

            cprint("\t Handling T=%.0fK, %d"%(T,len(props_tarr)))
            for pi,props in enumerate(props_tarr):
                i,phi = first(phi_arr,lambda phi: phi==props['phi'])
                p,P = first(P_arr,lambda P: P==props['P'])
                AX[i,p].bar(x-0.2,tdata_tarr[pi][rank],width=0.4)
                AX[i,p].bar(x+0.2,sdata_tarr[pi][rank],width=0.4)
                # if i==2 and p==1: plt.xticks(x,xtick,rotation=45)
                err_arr.append(1-np.dot(tdata_tarr[pi],sdata_tarr[pi]))
            set_sub_plots(AX, xlabel=r'sensitive reactions', ylabel=r'sensitivity',
                    legend=["brute force","adjoint"],ylim=[-1.,1.])
            save_figure(fig, path=figs_dir+'/sens_comp/%s_T=%.0fK.eps'%(name,T))
        cprint("Mean err: %.3e, Max err: %.3e"%(np.mean(err_arr), np.max(err_arr)), 'g')
    plt.show()
