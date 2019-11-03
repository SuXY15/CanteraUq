from utils import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

name = "H2"
mech = "mech/H2/H2.cti"
# phi_arr=[0.5]
# P_arr = [10.]
# T_arr = [650.]
# conditions = get_conditions(phi_arr, P_arr, T_arr)

props['mri'] = []
# calculator
def calculator(name, gas, props_arr, file_name):
    data_dict = {'props':[], 'tdata':[]}
    NUM = len(props_arr)
    for i,props in enumerate(props_arr):
        print("|{0:07.3f}s| calculator {1}: {2}/{3}...".format(time.time()-t0, name, i, NUM))
        # set all multipliers
        gas.set_multiplier(1.0) 
        for j,factor in enumerate(props['factor']):
            gas.set_multiplier(factor, j)

        # brute force method for idt sensitivity
        sens = np.zeros(gas.n_reactions)
        idt = get_ign(gas, props)

        for p in props['mri']:
            gas.set_multiplier(props['factor'][p]*(1+pdiff), p)
            pidt = get_ign(gas, props)
            sens[p] = (np.log(pidt) - np.log(idt))/np.log(1+pdiff)
            gas.set_multiplier(props['factor'][p], p)

        props['idt'] = idt
        # save data
        if (np.sum(np.isinf(sens))+np.sum(np.isnan(sens)))==0:
            data_dict['props'] = deepcopy(props)
            data_dict['tdata'] = deepcopy(sens.tolist())
            acquireLock(maxt, size, rank)
            json_writer(file_name, data_dict)
        else:
            cprint("NaN or Inf occured at %s."%file_name)

def sampling():
    # load mech and uncertainty
    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech[:-3]+'txt', UF=UF)
    
    # load main reactions
    main_reaction_idx = set(json_reader(sens_dir+name+"_mr.json")[0]['mr'])
    props['mri'] = [int(i) for i in main_reaction_idx]
    for props['phi'],props['P'],props['T'] in conditions:
        file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                        UF,props['phi'],props['P'],props['T'],samplingSize)
        if checkexists(file_name, delete=False):
            cprint("File %s exists."%file_name, 'b')

        # sample for factors
        props_arr = []
        for s in range(samplingSize):#  monte carlo
            props['factor'] = list(get_factors(uncetainty_factors))
            props_arr.append(deepcopy(props))

        # calculations
        print("\nStart c%d of phi=%.1f P=%.1fatm T=%.1fK"%(rank,props['phi'],props['P'],props['T']))
        l,i = [int(li) for li in np.linspace(0, len(props_arr), num=size+1)],rank
        calculator('c%d'%i, gas, props_arr[l[i]:l[i+1]], file_name)

def show(show_flag=0):
    m,name,mech = mechs[0]
    props['phi'],props['P'],props['T'] = 1.0, 10.0, 1000.0

    # loading
    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]
    
    idx = [i for i,td in enumerate(tdata_arr) if np.sum(np.isinf(td)+np.isnan(td))==0]
    tdata_arr = [tdata_arr[i] for i in idx]
    props_arr = [props_arr[i] for i in idx]
    idt = np.array([props['idt'] for props in props_arr])
    print(file_name)
    print("Num of samples:",len(sens_data_arr), "Useful:", len(idx))

    # calculating C matrix
    tdata_arr = normalize(tdata_arr)
    C = np.transpose(tdata_arr) @ np.array(tdata_arr)

    # SVD decomposing
    U,S,VT = np.linalg.svd(C)

    # dimensional projection
    w1x = [np.dot(VT[0], 3.*np.log(props['factor'])/np.log(uncetainty_factors)) for props in props_arr]
    w2x = [np.dot(VT[1], 3.*np.log(props['factor'])/np.log(uncetainty_factors)) for props in props_arr]

    # figure(1): eigenvalue
    graph_rank = np.arange(1,21)
    tick = np.arange(0,21,2)
    fig1 = plt.figure("Eigenvalue",figsize=c2i(12,9))
    plt.plot(graph_rank,S[graph_rank-1], marker='o', markerfacecolor='none', color='k')
    plt.xticks(tick)
    plt.yscale(r'log')
    plt.xlabel(r'Eigenvalue Number')
    plt.ylabel(r'Eigenvalue')
    save_figure(fig1, path=figs_dir+'acts_Eigenvalue.png')
    
    # figure(2): one dimensional projection
    fig2 = plt.figure(r"$\mathbf{w}_1$ projection",figsize=c2i(12,9))
    plt.scatter(w1x, np.log10(idt), marker='o', color='', edgecolors='k')
    plt.xlabel(r'$\mathbf{w}_1^T \mathbf{x}$')
    plt.ylabel(r'$\log_{10}(IDT)$')
    save_figure(fig2, path=figs_dir+'acts_W1projection.png')
    
    # figure(3): w1 components
    x,y = np.transpose([[i,v] for i,v in enumerate(VT[0]) if abs(v)>0.05])
    
    eqs = [r.equation for r in gas.reactions()]
    sVT = sorted([(i,v) for i,v in enumerate(VT[0])], key=lambda vi:-abs(vi[1]))

    fig3 = plt.figure("w1 components",figsize=c2i(12,9))
    markerline, stemlines, baseline = plt.stem(x,y, markerfmt='ko', linefmt='k-.',basefmt='gray')
    plt.plot([1,gas.n_reactions],[0,0],'-',color='gray')
    plt.setp(markerline, color='k', markerfacecolor='none', linewidth=2)
    plt.xlabel(r'Reaction Index')
    plt.ylabel(r'$\mathbf{w}_1$ components')
    plt.xlim([1, gas.n_reactions])
    
    for i in range(3):
        print(i, sVT[i][0], sVT[i][1], eqs[sVT[i][0]])
    
    save_figure(fig3, path=figs_dir+'acts_W1components.png')
    
    # figure(4): histogram
    # fig4 = plt.figure("histogram",figsize=c2i(12,9))
    # plt.title("histogram of IDT distribution")
    fig4, ax = plt.subplots(figsize=c2i(12,9))
    # plt.hist(np.log10(idt),bins=20,histtype='step')
    import pandas as pd
    dist = pd.DataFrame(np.log10(idt))
    dist.plot.kde(ax=ax, legend=False, color='k', lw=1)
    dist.plot.hist(density=True, ax=ax, bins=32, color='k', histtype='step', legend=False)

    plt.xlabel(r'$\log_{10}({IDT}[s])$')
    plt.ylabel(r'Normalized Histogram')
    save_figure(fig4, path=figs_dir+'acts_Histogram.png')
    
    # # figure(5): scatter
    # fig5 = plt.figure("w1 and w2 projection",figsize=c2i(12,9))
    # plt.scatter(w1x, w2x, c=np.log10(idt), s=50.0, cmap='rainbow',
    #             vmin=np.min(np.log10(idt)), vmax=np.max(np.log10(idt)))
    # plt.title("w1 and w2 projection of subspace")
    # plt.xlabel(r'$w_1^TX$')
    # plt.ylabel(r'$w_2^TX$')
    # save_figure(fig5, path=figs_dir+'acts_scatter.png')

    plt.show()

def generate():
    for m,name,mech in mechs:
        # loading
        cprint("Generating %s"%name, 'g')
        gas = load_mech(mech)
        uncetainty_factors = load_uncertainty(mech, UF=UF)
        for i,condition in enumerate(conditions):
            if i%size != rank: continue
            
            props['phi'],props['P'],props['T'] = condition
            file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                    UF,props['phi'],props['P'],props['T'],samplingSize)
            subs_name = file_name[:-5]+"_as.json"

            if not checkexists(subs_name,delete=False):
                print("Generating %s."%subs_name)
                sens_data_arr = json_reader(file_name)
                print("Num of samples:",len(sens_data_arr))
                tdata_arr = [sd['tdata'] for sd in sens_data_arr]
                props_arr = [sd['props'] for sd in sens_data_arr]
                idt = np.array([props['idt'] for props in props_arr])

                idx = [i for i,td in enumerate(tdata_arr) if np.sum(np.isinf(td)+np.isnan(td))==0]
                tdata_arr = [tdata_arr[i] for i in idx]
                props_arr = [props_arr[i] for i in idx]
                # calculating C matrix
                
                tdata_arr = normalize(tdata_arr)
                C = np.transpose(tdata_arr) @ np.array(tdata_arr)
                
                try:
                    # SVD decomposing
                    U,S,VT = np.linalg.svd(C)
                    cprint("S[0]:%.4f,S[1]:%.4f,S[2]:%.4f"%(S[0],S[1],S[2]),'b')
                    json_writer(subs_name,{'subspace':VT.tolist()})
                except:
                    cprint("Error in SVD decomposing")

# main
if __name__=="__main__":
    checkpath(acts_dir)
    if len(sys.argv) < 3:
        cprint("Input save or show in args!")
        exit(0)
    if sys.argv[2] == 'sampling':
        sampling()
    if sys.argv[2] == 'show':
        show()
    if sys.argv[2] == 'generate':
        generate()