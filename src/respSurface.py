from utils import *
sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *

def fitting_plots(x, xf, xp, v, vf, vp, NX, NP):
    fig1 = plt.figure("Training and Validation",figsize=c2i(12,6))
    # plt.title("Training and Validation")
    plt.plot(xf,xp,'kv',ms=6,fillstyle='none')
    plt.plot(vf,vp,'r^',ms=6,fillstyle='none')
    fmin, fmax = np.min(xf), np.max(xf)
    plt.plot([fmin,fmax],[fmin,fmax],'k--')
    plt.xlabel(r'$\log_{10}(\tau_{train})$')
    plt.ylabel(r'$\log_{10}(\tau_{prediction})$')
    plt.legend(['train','valid'], loc='lower right',bbox_to_anchor=(1.0,-0.1),frameon=False)

    fig2 = plt.figure("Response surface prediction",figsize=c2i(12,6))
    # plt.title("Response surface prediction")
    plt.plot(NX[:,0],NP,'k.',ms=.4)
    plt.plot(x[:,0],xp,'kv',ms=6,fillstyle='none')
    plt.plot(v[:,0],vp,'r^',ms=6,fillstyle='none')

    plt.xlabel(r'$\mathbf{w}_1^T \mathbf{x}$')
    plt.ylabel(r'$\log_{10}(\rm{IDT}[s])$')
    plt.ylim([-4,-1])
    plt.legend(["samples","train","valid"], loc='lower right', bbox_to_anchor=(1.0,-0.1),frameon=False)
    return fig1, fig2

def load_training_set(file_name):
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]

    # prepare data
    idx = [i for i,props in enumerate(props_arr) if props['idt']>0]
    idt = np.array([props_arr[i]['idt'] for i in idx])
    tdata = np.array([tdata_arr[i] for i in idx])
    factor = np.array([props_arr[i]['factor'] for i in idx])
    return idt, factor, tdata

def compare_gradients(pdx, rdx):
    rdx = rdx/np.max(rdx)
    pdx = pdx/np.max(pdx)
    cos = [cosine(rdx[i],pdx[i]) for i in range(len(rdx))]
    co2 = [cosine(rdx[:,d],pdx[:,d]) for d in range(len(rdx[0]))]
    print(pdx.shape,rdx.shape, "%.6f"%(co2[0]))

    fig = plt.figure("Gradients Comparison")
    plt.title("Gradients Comparison")
    for i in range(len(rdx[0])):
        plt.subplot(121);plt.plot(pdx[:,i],rdx[:,i],'%s'%symbol_arr[i],fillstyle='none')
        plt.subplot(122);
        plt.plot(pdx[:,i],'b%s'%symbol_arr[i],fillstyle='none')
        plt.plot(rdx[:,i],'r%s'%symbol_arr[i],fillstyle='none')
    plt.plot(cos,'C9')
    return fig

def train(dim=3,order=3,N=50000,method=''):
    m,name,mech = mechs[0]
    props = config_dict['props']
    props['phi'], props['P'], props['T'] = 1.0, 1.0, 1000.
    if method == 'ann':
        response_surface = NerualNetworkApproximation(dim=dim, N=10)
    else:
        response_surface = PolynomialApproximation(N=order)
        # response_surface = RadialBasisApproximation(N=order)

    # prepare data
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
    idt, factor, tdata = load_training_set(file_name)
    
    # # training
    # f, X = [], []
    # for j in range(len(idt)):
    #     iidt, fact, sens = idt[j], factor[j], tdata[j]
    #     for i in range(len(fact)):
    #         if not np.isinf(sens[i]) and not np.isnan(sens[i]) and sens[i]!=0:
    #             fact[i] *= (1+pdiff)
    #             f.append(np.log10(np.exp(sens[i]*np.log(1+pdiff)+np.log(iidt))))
    #             X.append(np.dot(VT[:dim], 3.*np.log(fact)/np.log(uncetainty_factors)))
    #             fact[i] /= (1+pdiff)
    # X = np.array(X)
    # f = np.array(f).reshape(len(f),1)
    # trainSize = len(f)>>1
    # x,v = X[:trainSize,:],X[trainSize:,:]
    # xf,vf = f[:trainSize,:],f[trainSize:,:]
    # if method=='ann':
    #     response_surface.train(x,xf,v,vf)
    # else:
    #     response_surface.train(x,xf)
    
    f = np.transpose([np.log10(idt),])
    X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
    trainSize = len(f)>>1
    x,v = X[:trainSize,:],X[trainSize:,:]
    xf,vf = f[:trainSize,:],f[trainSize:,:]
    if method=='ann':
        response_surface.train(x,xf,v,vf)
    else:
        response_surface.train(x,xf)

    # predicting
    xp,dx = response_surface.predict(x, compgrad=True)
    vp,dv = response_surface.predict(v, compgrad=True)

    # check gradients
    pdx = np.vstack((dx,dv))
    rdx = np.dot(tdata,np.transpose(VT[:dim]))
    compare_gradients(pdx,rdx)

    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
    NP,dn = response_surface.predict(NX)
    
    # showing
    print(file_name)
    cprint("dim: %d, ord: %d, train: %d, valid: %d"%(dim,order,trainSize,len(f)-trainSize), 'b')
    if method=='ann': response_surface.show_loss()
    fig1,fig2 = fitting_plots(x, xf, xp, v, vf, vp, NX, NP)
    
    import pandas as pd
    fig3, ax = plt.subplots(figsize=c2i(12,9))

    dist = pd.DataFrame(NP, columns=['samples'])
    dist.plot.kde(ax=ax, legend=False, color='k', lw=1)
    dist.plot.hist(density=True, ax=ax, bins=32, color='k', histtype='step')

    dist = pd.DataFrame(xf.flatten().tolist()+vf.flatten().tolist(), columns=['train'])
    dist.plot.kde(ax=ax, legend=False, color='r', lw=1)
    dist.plot.hist(density=True, ax=ax, bins=32, color='r', histtype='step')
    
    dist = pd.DataFrame(xp.flatten().tolist()+vp.flatten().tolist(), columns=['predict'])
    dist.plot.kde(ax=ax, legend=False, color='b', lw=1)
    dist.plot.hist(density=True, ax=ax, bins=32, color='b', histtype='step')

    plt.xlim([-2.3, -0.6])
    plt.xlabel(r'$\log_{10}(\tau)$')
    plt.ylabel(r'Normalized Histogram')
    plt.legend(["samples", "train", "predict"], frameon=False)
    save_figure(fig1, figs_dir+"resp_train&valid.eps")
    save_figure(fig2, figs_dir+"resp_prediction.pdf")
    save_figure(fig3, figs_dir+"resp_histogram.eps")
    plt.show()

# train without subspace
def rawtrain(order=2,N=50000):
    m,name,mech = mechs[0]
    props = config_dict['props']
    props['phi'], props['T'], props['P'] = 0.5, 1600., 1.0
    response_surface = PolynomialApproximation(N=order)
    # response_surface = RadialBasisApproximation(N=order)

    # loading file
    gas = load_mech(mech)
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    sens_data_arr = json_reader(file_name)
    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
    props_arr = [sd['props'] for sd in sens_data_arr]
    mri = props_arr[0]['mri']
    dim = len(mri)

    idt = np.array([props['idt'] for props in props_arr])
    idx = [i for i,iidt in enumerate(idt) if iidt>0]
    # f = np.transpose([np.log10(idt[idx]),])
    # X = np.array([(3.*np.log(p['factor'])/np.log(uncetainty_factors))[mri] for p in props_arr])[idx]
    f, X = [], []
    for sd in sens_data_arr:
        props, tdata = sd['props'],sd['tdata']
        factor = props['factor']
        for i in props['mri']:
            if not np.isinf(tdata[i]) and not np.isnan(tdata[i]) and props['idt']>0:
                factor[i] *= (1+pdiff)
                f.append(np.log10(np.exp(tdata[i]*np.log(1+pdiff)+np.log(props['idt']))))
                X.append((3.*np.log(factor)/np.log(uncetainty_factors))[mri])
                factor[i] /= (1+pdiff)
    X = np.array(X)
    f = np.array(f).reshape(len(f),1)

    # training
    trainSize = len(f)>>1
    x,v = X[:trainSize,:],X[trainSize:,:]
    xf,vf = f[:trainSize,:],f[trainSize:,:]
    response_surface.train(x,xf)

    # predicting
    xp,dx = response_surface.predict(x)
    vp,dv = response_surface.predict(v)
    NX = np.transpose(np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))[mri])
    NP,dn = response_surface.predict(NX)

    # showing
    print(file_name)
    cprint("dim: %d, ord: %d, train: %d, valid: %d"%(dim,order,trainSize,len(f)-trainSize), 'b')
    fig1,fig2 = fitting_plots(x, xf, xp, v, vf, vp, NX, NP)
    plt.show()

def distribution(dim=3,N=50000):
    for m,name,mech in mechs:
        cprint("loading mech: %s"%mech, 'g')
        uncetainty_factors = load_uncertainty(mech, UF=UF)
        # Load data
        data_dict, dist_dict = {}, {}
        maxA ,minA = -1e10, 1e10
        for idx_phi,phi in enumerate(phi_arr):
            data_dict[phi], dist_dict[phi] = {}, {}
            
            for idx_P,P in enumerate(P_arr):
                data_list, dist_list = [], []
                
                for idx_T,T in enumerate(T_arr):
                    props['phi'],props['P'],props['T'] = phi,P,T
                    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                        UF,props['phi'],props['P'],props['T'],samplingSize)
                    sens_data_arr = json_reader(file_name)
                    tdata_arr = [sd['tdata'] for sd in sens_data_arr]
                    props_arr = [sd['props'] for sd in sens_data_arr]
                    subs_name = file_name[:-5]+"_as.json"
                    # prepare training
                    VT = np.array(json_reader(subs_name)[0]['subspace'])
                    idt = np.array([props['idt'] for props in props_arr])
                    idx = [i for i,iidt in enumerate(idt) if iidt>0]
                    f = np.transpose([np.log10(idt[idx]),])
                    X = np.array([np.dot(VT[:dim],3.*np.log(props_arr[i]['factor'])\
                                /np.log(uncetainty_factors)) for i in idx])

                    print("Data used:",len(f),"throwed:",len(idt)-len(idx))
                    minA, maxA = min(minA, np.mean(f)-1), max(maxA, np.mean(f)+1)
                    data_list.append(np.mean(f))
                    dist_list.append([X[:,0],f])
                data_dict[phi][P] = data_list
                dist_dict[phi][P] = dist_list

        # prepare figures
        figA, AX = get_sub_plots(num="Data Distribution, %s"%name)
        minT, maxT = 1000./np.max(T_arr), 1000./np.min(T_arr)
        minA, maxA = maxA-(maxA-minA)*6/5, (maxA-minA)*6/5+minA
        minT, maxT = maxT-(maxT-minT)*19/18, (maxT-minT)*19/18+minT
        T_revert = np.array([1000./Ti for Ti in T_arr])
        for i,phi in enumerate(phi_arr):
            for p,P in enumerate(P_arr):
                data_list = data_dict[phi][P]
                dist_list = dist_dict[phi][P]
                AX[i][p].plot(T_revert, data_list,'k-')
                for t,T in enumerate(T_arr):
                    height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                    b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                    width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                    l_width = (0-min(dist_list[t][0]))/6*0.03
                    left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                    bottom = (2-i)*0.3 + 0.05 + (data_list[t]-minA)/(maxA-minA)*0.3-b_height
                    ax = plt.axes([left,bottom,width,height],facecolor='none')
                    plt.plot(dist_list[t][0],dist_list[t][1],'ko',ms=1,fillstyle='full')
                    for edge in ['top','right','bottom','left']:
                        ax.spines[edge].set_visible(False)
                    plt.xticks([])
                    plt.yticks([])

        set_sub_plots(AX, r'1000/T, $K^{-1}$', r'$\log_{10}(\tau)$', [name],
                            xlim=[minT,maxT],ylim=[minA,maxA])
        save_figure(figA, path=figs_dir+'resp_Distribution_%s.pdf'%name)
    plt.show()

def predict(dim=3,order=2,N=50000):
    # run calculation
    response_surface = PolynomialApproximation(N=order)
    # response_surface = RadialBasisApproximation(N=order)

    for m,name,mech in mechs:
        gas = load_mech(mech)
        resp_name = resp_dir+"%s_dim=%d_order=%d.json"%(name,dim,order)
        if checkexists(resp_name,delete=False):
            print("File %s exists."%resp_name)
        else:
            data_dict = {}
            for idx_phi,phi in enumerate(phi_arr):
                data_dict[phi] = {}
                for idx_P,P in enumerate(P_arr):
                    data_list = []
                    for idx_T,T in enumerate(T_arr):
                        # preparation
                        props['phi'],props['P'],props['T'] = phi,P,T
                        file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                                    UF,props['phi'],props['P'],props['T'],samplingSize)
                        uncetainty_factors = load_uncertainty(mech, UF=UF)
                        print("Running training and prediction for %s"%file_name)

                        # loading data
                        sens_data_arr = json_reader(file_name)
                        tdata_arr = [sd['tdata'] for sd in sens_data_arr]
                        props_arr = [sd['props'] for sd in sens_data_arr]
                        subs_name = file_name[:-5]+"_as.json"
                        VT = np.array(json_reader(subs_name)[0]['subspace'])

                        # training
                        idt = np.array([props['idt'] for props in props_arr])
                        idx = [i for i,iidt in enumerate(idt) if iidt>0]
                        f = np.transpose([np.log10(idt[idx]),])
                        X = np.array([np.dot(VT[:dim],3.*np.log(props_arr[i]['factor'])\
                                    /np.log(uncetainty_factors)) for i in idx])

                        print("Data used:",len(f),"throwed:",len(idt)-len(idx))
                        try:
                            response_surface.train(X,f)
                        except:
                            print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                            data_list.append([None,None,None])
                            continue
                        # predicting
                        NX = np.transpose(np.dot(VT[:dim],np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
                        NP,dn = response_surface.predict(NX)
                        data_list.append([np.mean(NP),np.std(NP),len(sens_data_arr)])
                    data_dict[phi][P] = data_list
            json_writer(resp_name,data_dict)

    # prepare figures
    maxA, minA = -1e10, 1e10
    maxB, minB = -1e10, 1e10
    maxC, minC = -1e10, 1e10
    figA, AX = get_sub_plots(num="IDT Error Comparison")
    figB, BX = get_sub_plots(num="Sigma Comparison")
    figC, CX = get_sub_plots(num="Sample Check")
    
    # m,name,mech = mechs[3]
    # gas = load_mech(mech)
    # resp_name = resp_dir+"%s_dim=%d_order=%d.json"%(name,dim,order)
    # data_dict = json_reader(resp_name)[0]
    # for p,P in enumerate(P_arr):
    #     DATA_list = []
    #     for i,phi in enumerate(phi_arr):
    #         data_list = []
    #         for data in data_dict[str(phi)][str(P)]:
    #             if data[0]!=None:
    #                 data_list.append(data[1])
    #             else:
    #                 data_list.append(np.nan)
    #         DATA_list.append(data_list)
    #     print("P:",P)
    #     print(np.array(DATA_list).T)

    # plotting
    for m,name,mech in mechs:
        cprint("Plotting %s"%name,'g')
        gas = load_mech(mech)
        resp_name = resp_dir+"%s_dim=%d_order=%d.json"%(name,dim,order)
        data_dict = json_reader(resp_name)[0]
        for i,phi in enumerate(phi_arr):
            for p,P in enumerate(P_arr):
                data_list = []
                for data in data_dict[str(phi)][str(P)]:
                    if data[0]!=None:
                        data_list.append(data)
                        maxA,minA = max(data[0]+2*data[1],maxA),min(data[0]-2*data[1],minA)
                        maxB,minB = max(data[1]*1.2,maxB),min(data[1]*0.8,minB)
                        maxC,minC = max(data[2]*1.2,maxC),min(data[2]*0.8,minC)
                    else:
                        data_list.append([np.nan,np.nan,np.nan])
                data_list = np.array(data_list)
                AX[i][p].errorbar(1000./np.array(T_arr),data_list[:,0],yerr=data_list[:,1],
                    fmt=color_arr[m]+"-"+symbol_arr[m])
                BX[i][p].plot(1000./np.array(T_arr),data_list[:,1],color_arr[m]+"-"+symbol_arr[m])
                CX[i][p].plot(1000./np.array(T_arr),data_list[:,2],color_arr[m]+"-"+symbol_arr[m])
    # figure setting
    set_sub_plots(AX, r'1000/T (K$^{-1}$)', r'$\log_{10}(\tau)$', mech_arr, ylim=[minA,maxA])
    set_sub_plots(BX, r'1000/T (K$^{-1}$)', r'$\sigma_r$', mech_arr, ylim=[minB,maxB])
    set_sub_plots(CX, r'1000/T (K$^{-1}$)', r'$N_{train}$', mech_arr, ylim=[minC,maxC])
    save_figure(figA, figs_dir+'resp_IDT.eps')
    save_figure(figB, figs_dir+'resp_sigma.eps')
    save_figure(figC, figs_dir+'resp_check.eps')
    plt.show()

def single_condition_propagation(dim=3,order=2,N=50000):
    ns_arr = []
    mu_arr = []
    sigma_arr = []
    rt_arr = []
    peqs = None
    seqs = None
    pWd = None
    sWd = None
    pRS = None
    rt = None
    for m,name,mech in mechs:
        props['phi'], props['P'], props['T'] = 1.0, 1.0, 1200.
        response_surface = PolynomialApproximation(N=order)

        # prepare data
        file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                    UF,props['phi'],props['P'],props['T'],samplingSize)
        uncetainty_factors = load_uncertainty(mech, UF=UF)
        VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
        Wd = np.transpose(VT[:dim])

        idt, factor, tdata = load_training_set(file_name)
        
        f = np.transpose([np.log10(idt),])
        X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
        trainSize = len(f)
        x,v = X[:trainSize,:],X[trainSize:,:]
        xf,vf = f[:trainSize,:],f[trainSize:,:]
        response_surface.train(x, xf)

        # # predicting
        # xp,dx = response_surface.predict(x, compgrad=True)
        # vp,dv = response_surface.predict(v, compgrad=True)

        # # check gradients
        # pdx = np.vstack((dx,dv))
        # rdx = np.dot(tdata,np.transpose(VT[:dim]))

        # NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
        # NP,dn = response_surface.predict(NX)
        
        # mu, sigma = np.mean(NP), np.std(NP)

        mu, sigma = response_surface.mean(Wd), response_surface.std(Wd)[0][0]
        
        gas = ct.Solution(mech)
        if peqs == None:
            peqs = [r.equation for r in gas.reactions()]
            pWd = Wd
            pRS = deepcopy(response_surface)
            rt_arr.append(None)
        else:
            seqs = [r.equation for r in gas.reactions()]
            IDR = parentRank(range(len(seqs)),seqs,peqs)
            S = np.zeros([len(peqs), len(seqs)]) # d x r
            for i,idr in enumerate(IDR):
                if(i>0 and IDR[i-1]==idr):
                    S[idr+1,i] = 1
                else:
                    S[idr,i] = 1
            Wi = S.transpose() @ pWd
            sigma_d = pRS.std(pWd)[0][0]
            sigma_t = pRS.std(S @ Wi)[0][0]
            rt = abs(sigma_d-sigma_t)/(abs(sigma_d-sigma_t)+abs(sigma_t-sigma))
            pRS = deepcopy(response_surface)
            pWd = deepcopy(Wd)
            peqs = deepcopy(seqs)
            rt_arr.append(rt)

        ns_arr.append(gas.n_species)
        mu_arr.append(mu)
        sigma_arr.append(sigma)
        print(mu, sigma, rt)
    fig = plt.figure(figsize=c2i(12,9))
    ax = plt.subplot(111)
    ax.plot(ns_arr, mu_arr, color_arr[1]+line_arr[1]+symbol_arr[1], fillstyle='none')
    ax2 = ax.twinx()
    ax2.plot(ns_arr, sigma_arr, color_arr[2]+line_arr[2]+symbol_arr[2], fillstyle='none')
    ax2.plot(ns_arr, rt_arr, 'k-')
    ax.set_xlabel("Number of Species")
    ax.set_ylabel(r"$\mu_r$")
    ax2.set_ylabel(r"$\sigma_r$")
    save_figure(fig, path=figs_dir+"single_condtion.png")
    plt.show()

if __name__ == '__main__':
    checkpath(resp_dir)
    if len(sys.argv) < 3:
        print("Input train or predict in args!")
        exit(0)
    if sys.argv[2] == 'predict':
        predict(dim=3,order=2,N=50000)
    if sys.argv[2] == 'distribution':
        distribution()
    if sys.argv[2] == 'train':
        train()
    if sys.argv[2] == 'rawtrain':
        train(method='raw')
    if sys.argv[2] == 'anntrain':
        train(method='ann')
    if sys.argv[2] == 'single':
        single_condition_propagation()