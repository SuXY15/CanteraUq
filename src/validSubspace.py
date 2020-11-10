from utils import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

name = "DME55"
mech = "mech/DME/DMEzhao.cti"

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

samplingSize=6400

def sampling():
    # load mech and uncertainty
    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    props['mri'] = [int(i) for i in range(gas.n_reactions)]

    props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200
    file_name = acts_dir+"valid_%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                    UF,props['phi'],props['P'],props['T'],samplingSize)
    if checkexists(file_name, delete=False):
        cprint("File %s exists."%file_name, 'b')

    # sample for factors
    props_arr = []
    for s in range(400):#  monte carlo
        props['factor'] = list(get_factors(uncetainty_factors))
    for s in range(1200):#  monte carlo
        props['factor'] = list(get_factors(uncetainty_factors))
        props_arr.append(deepcopy(props))

    # calculations
    print("\nStart c%d of phi=%.1f P=%.1fatm T=%.1fK"%(rank,props['phi'],props['P'],props['T']))
    l,i = [int(li) for li in np.linspace(0, len(props_arr), num=size+1)],rank
    calculator('c%d'%i, gas, props_arr[l[i]:l[i+1]], file_name)

sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *
response_surface = PolynomialApproximation(N=2)

def show(show_flag=0):
    props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200.0
    dim, N = 3, 50000

    # get sensitivity
    sens_name = sens_dir+"DMEzhao.json"
    for data in json_reader(sens_name):
        dp = data['props']
        if dp['phi']==props['phi'] and dp['T']==props['T'] and dp['P']==props['P']:
            sens_data = normalize(data['tdata'])
    ridx = np.argsort(np.abs(sens_data))

    # loading data
    gas = load_mech(mech)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    file_name = acts_dir+"valid_%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]
    
    idx = [i for i,td in enumerate(tdata_arr) if np.sum(np.isinf(td)+np.isnan(td))==0]
    tdata = np.array([tdata_arr[i] for i in idx])
    props_arr = np.array([props_arr[i] for i in idx])
    idt = np.array([props_arr[i]['idt'] for i in idx]) # mx1
    factors = np.array([3.*np.log(props_arr[i]['factor'])/np.log(uncetainty_factors) for i in idx]) # mxn
    tDATA = deepcopy(tdata)

    # calculating subspace
    def cal_subspace(tdata):       
        tdata = normalize(tdata)
        C = np.transpose(tdata) @ np.array(tdata)
        U,S,VT = np.linalg.svd(C)
        return VT

    def cal_responseSurface(idt, factors, VT):
        # train response surface
        f = np.transpose([np.log10(idt),])
        X = np.array([np.dot(VT[:dim],factor) for factor in factors])
        response_surface.train(X,f)

        # predicting
        # NX = np.transpose(np.dot(VT[:dim],np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
        # NP,dn = response_surface.predict(NX)
        # return np.mean(NP),np.std(NP)

        Wd = np.transpose(VT[:dim])
        mean_v, std_v = response_surface.mean(Wd), response_surface.std(Wd)
        return mean_v, std_v[0][0]

    td0 = deepcopy(tdata)
    VT0 = cal_subspace(tdata)
    mu0, sigma0 = cal_responseSurface(idt, factors, VT0)
    td0norm = [np.dot(td0[i], td0[i]) for i in range(len(td0))]

    sens_ratio_arr = []
    subs_ratio_arr = []
    mu_ratio_arr = []
    sigma_ratio_arr = []
    rawsigma_ratio_arr = []

    import random
    oidx = range(len(idt))
    for k in range(35, 285):
        # progress(k/100)
        # calculating subspace
        tdata[:,ridx[:k]] = 0 # set the most in-sensitive to be zero
        sens_ratio = np.mean([np.dot(td0[i], tdata[i])/td0norm[i] for i in range(len(td0))])

        if len(sens_ratio_arr)>0 and (1-sens_ratio)/(1-sens_ratio_arr[-1]) < 1.02:
            continue

        for z in range(10):
            sidx = random.sample(oidx, 400)
            # sidx = oidx
            VT = cal_subspace(tdata[sidx])
            mu, sigma = cal_responseSurface(idt[sidx], factors[sidx], VT)
            
            sens_ratio_arr.append(sens_ratio)
            subs_ratio_arr.append(np.dot(VT[0],VT0[0]))
            mu_ratio_arr.append(mu/mu0)
            sigma_ratio_arr.append(sigma/sigma0)

            # rawsigma = np.std(np.log10(idt[sidx]))
            # rawsigma_ratio_arr.append(rawsigma/sigma0)

            if sens_ratio>0.99 and abs(1-sigma/sigma0)>2e-2:
                print("%.5f %.3e"%(sens_ratio, abs(1-sigma/sigma0)))

    sens_ratio_arr = 1-np.array(sens_ratio_arr)
    subs_ratio_arr = 1-np.array(subs_ratio_arr)
    mu_ratio_arr = np.abs(1-np.array(mu_ratio_arr))
    sigma_ratio_arr = np.abs(1-np.array(sigma_ratio_arr))
    # rawsigma_ratio_arr = np.abs(1-np.array(rawsigma_ratio_arr))

    fig = plt.figure(figsize=c2i(12,12))
    ax2 = plt.subplot(212)
    ax = plt.subplot(211)
    ax.plot(sens_ratio_arr, subs_ratio_arr, color_arr[0]+symbol_arr[0], fillstyle='none',ms=4, alpha=0.5)
    ax.plot(sens_ratio_arr, mu_ratio_arr, color_arr[1]+symbol_arr[1], fillstyle='none',ms=4, alpha=0.5)
    ax.plot(sens_ratio_arr, sigma_ratio_arr, color_arr[2]+symbol_arr[2], fillstyle='none',ms=4, alpha=0.5)
    # ax.plot(sens_ratio_arr, rawsigma_ratio_arr, color_arr[3]+symbol_arr[3], fillstyle='none',ms=4, alpha=0.5)
    ax.plot([np.min(sens_ratio_arr),np.max(sens_ratio_arr)], [2.5e-2, 2.5e-2], 'k--')
    ax.set_ylim([3e-6, 2e-1])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel("Prediction Error")
    plt.legend([r"$w_1$", r'$\mu_{r}$', r'$\sigma_{r}$'], loc=(0.7, -0.8))

    filename = sens_dir+"DMEzhao"+".json"
    sens_data = [d['tdata'] for d in json_reader(filename)]
    filename = sens_dir+"DMEzhao"+"_mr.json"
    mri = json_reader(filename)[0]['mr']
    
    # errors in sens data at nominal parameters, for all conditions
    sens_error_arr = [] 
    for sens in sens_data:
        orig_s = normalize(sens)
        main_s = np.zeros(len(orig_s))
        main_s[mri] = orig_s[mri]
        ratio = np.dot(orig_s, main_s)
        sens_error_arr.append(1-ratio)

    # errors in the parameter space, for the single condition
    sens_error_arr = []
    for sens in tDATA:
        orig_s = normalize(sens)
        main_s = np.zeros(len(orig_s))
        main_s[mri] = orig_s[mri]
        ratio = np.dot(orig_s, main_s)
        sens_error_arr.append(1-ratio)

    ax2.hist(sens_error_arr, alpha=0.5, bins=np.logspace(np.log10(np.min(sens_ratio_arr)),np.log10(np.max(sens_ratio_arr)), 30))
    ax2.set_ylabel("PDF of Sensitvity Error")
    ax2.set_xscale('log')
    ax2.set_xlabel("Sensitivity Error")
    save_figure(fig, path=figs_dir+"validSubspace_show.png")
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