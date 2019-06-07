from utils import *
sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *

dim, order, N = 3,2,50000

Type = int(sys.argv[2])
Flag = None
Legend = []
TypeName = ['O','I', 'II', 'III', 'IV']

if __name__ == "__main__":
    pmech = mechs[1]
    smech = mechs[2]

    # preparing
    pgas = load_mech(pmech[2])
    sgas = load_mech(smech[2])
    peqs = [r.equation for r in pgas.reactions()]
    seqs = [r.equation for r in sgas.reactions()]
    p_uf = load_uncertainty(pmech[2], UF=UF)
    s_uf = load_uncertainty(smech[2], UF=UF)
    IDR = parentRank(range(len(seqs)),seqs,peqs)
    S = np.zeros([len(peqs), len(seqs)]) # d x r
    for i,idr in enumerate(IDR):
        if(i>0 and IDR[i-1]==idr):
            S[idr+1,i] = 1
        else:
            S[idr,i] = 1
    # print([[i for (i,si) in enumerate(Si) if si]for Si in S])
    # sys.exit()
    response_surface = PolynomialApproximation(N=order)

    # = = = = = = = = = = = =
    # pmech
    # Legend.append(pmech[1])
    m,name,mech = pmech
    cprint("Printing %s"%name,'g')

    # Load data
    data_dict = {}
    maxA ,minA = -1e10, 1e10
    maxB ,minB = -1e10, 1e10
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi] = {}

        for idx_P,P in enumerate(P_arr):
            data_list = []

            for idx_T,T in enumerate(T_arr):
                progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wd = get_subspace(pmech[1], props, dim)
                Xd, f = get_Xy(pmech[1], props, p_uf)
                X = Xd @ Wd
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wd

                try:
                    response_surface.train(X,f)
                    v,dv = response_surface.predict(NX)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                    print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                    data_list.append([None,None])
                    continue
            data_dict[phi][P] = data_list

    # prepare figures
    figB, BX = get_sub_plots(num="Data Sigmas")
    pmech_data = deepcopy(data_dict)

    Type = 3
    # = = = = = = = = = = = =
    # pmech -> intermediate
    Legend.append("transition - "+TypeName[Type])
    cprint("\nPrinting transition - "+TypeName[Type],'g')
    data_dict = {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi] = {}

        for idx_P,P in enumerate(P_arr):
            data_list = []

            for idx_T,T in enumerate(T_arr):
                progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wd = get_subspace(pmech[1], props, dim)
                Wi = normalize(S.transpose() @ Wd)
                Xd, f = get_Xy(pmech[1], props, p_uf)
                X = Xd @ S @ Wi
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi

                try:
                    response_surface.train(X,f)
                    v,dv = response_surface.predict(NX)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                   print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                   data_list.append([None,None])
                   continue
            data_dict[phi][P] = data_list
        
    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = (np.array(data_dict[phi][P]) - np.array(pmech_data[phi][P]))/np.array(pmech_data[phi][P])
            data_B = data_list[:,1]
            BX[i][p].plot(T_revert, data_B, color_arr[0])
            minB, maxB = min(minB, min(data_B)*1.2), max(maxB, max(data_B)*1.2)

    Type = 2
    # = = = = = = = = = = = =
    # pmech -> intermediate 2 if necessary
    Legend.append("transition - "+TypeName[Type])
    cprint("\nPrinting transition - "+TypeName[Type],'g')
    data_dict = {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi] = {}

        for idx_P,P in enumerate(P_arr):
            data_list = []

            for idx_T,T in enumerate(T_arr):
                progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wr = get_subspace(smech[1], props, dim)
                Wi = S @ Wr
                Xd, f = get_Xy(pmech[1], props, p_uf)
                X = Xd @ Wi
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wi

                try:
                    response_surface.train(X,f)
                    v,dv = response_surface.predict(NX)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                   print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                   data_list.append([None,None])
                   continue
            data_dict[phi][P] = data_list
            
    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = (np.array(data_dict[phi][P]) - np.array(pmech_data[phi][P]))/np.array(pmech_data[phi][P])
            data_B = data_list[:,1]
            BX[i][p].plot(T_revert, data_B, color_arr[1])
            minB, maxB = min(minB, min(data_B)*1.2), max(maxB, max(data_B)*1.2)

    # = = = = = = = = = = = =
    # smech
    Legend.append(smech[1])
    m,name,mech = smech
    cprint("\nPrinting %s"%name,'g')
    
    # Load data
    data_dict = {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi] = {}
        
        for idx_P,P in enumerate(P_arr):
            data_list = []

            for idx_T,T in enumerate(T_arr):
                progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T
                
                Wr = get_subspace(smech[1], props, dim)
                Xr, f = get_Xy(smech[1], props, s_uf)
                X = Xr @ Wr
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wr

                try:
                    response_surface.train(X,f)
                    v,dv = response_surface.predict(NX)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                    print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                    data_list.append([None,None])
                    continue
            data_dict[phi][P] = data_list

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = (np.array(data_dict[phi][P]) - np.array(pmech_data[phi][P]))/np.array(pmech_data[phi][P])
            data_B = data_list[:,1]
            BX[i][p].plot(T_revert, data_B, color_arr[2])
            minB, maxB = min(minB, min(data_B)*1.2), max(maxB, max(data_B)*1.2)
    
    minT, maxT = 1000./np.max(T_arr), 1000./np.min(T_arr)
    minT, maxT = maxT-(maxT-minT)*19/18, (maxT-minT)*19/18+minT
    set_sub_plots(BX, r'1000/T, $K^{-1}$', r'$\Delta \sigma/\sigma_p$',Legend, xlim=[minT,maxT],ylim=[minB,maxB])
    save_figure(figB, figs_dir+'prop_err_%s_%s_Sigm.png'%(pmech[1], smech[1]))
    plt.show()
