from utils import *
sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *

dim, order, N = 3,2,50000

Type = int(sys.argv[2])
Flag = None
Legend = []
TypeName = ['O','I', 'II', 'III', 'IV']

if __name__ == "__main__":
    pmech = mechs[2]
    smech = mechs[3]

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
    Legend.append(pmech[1])
    m,name,mech = pmech
    cprint("Printing %s"%name,'g')

    # Load data
    data_dict, dist_dict = {}, {}
    maxA ,minA = -1e10, 1e10
    maxB ,minB = -1e10, 1e10
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi], dist_dict[phi] = {}, {}
        
        for idx_P,P in enumerate(P_arr):
            data_list, dist_list = [], []

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
                    minA, maxA = min(minA, np.mean(v)), max(maxA, np.mean(v))
                    minB, maxB = min(minB, np.std(v))*0.8, max(maxB, np.std(v)*1.2)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                    print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                    data_list.append([None,None])
                    continue
                dist_list.append([X[:,0],f])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list

    # prepare figures
    figA, AX = get_sub_plots(num="Data Distributions")
    figB, BX = get_sub_plots(num="Data Sigmas")
    minT, maxT = 1000./np.max(T_arr), 1000./np.min(T_arr)
    minA, maxA = maxA-(maxA-minA)*6/5, (maxA-minA)*6/5+minA
    minT, maxT = maxT-(maxT-minT)*19/18, (maxT-minT)*19/18+minT
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            BX[i][p].plot(T_revert, data_list[:,1],'C0')
            AX[i][p].plot(T_revert, data_list[:,0],'C0')
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure(num="Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],'C0.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])

    if Type==5:
        Flag = True
        Type = 3
    # = = = = = = = = = = = =
    # pmech -> intermediate
    Legend.append("transition - "+TypeName[Type])
    cprint("\nPrinting transition - "+TypeName[Type],'g')
    data_dict, dist_dict = {}, {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi], dist_dict[phi] = {}, {}

        for idx_P,P in enumerate(P_arr):
            data_list, dist_list = [], []

            for idx_T,T in enumerate(T_arr):
                progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                if Type==1:
                    Wd = get_subspace(pmech[1], props, dim)
                    Wi = normalize(S.transpose() @ Wd)
                    Xr, f = get_Xy(smech[1], props, s_uf)
                    X = Xr @ Wi
                    NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
                if Type==2:
                    Wr = get_subspace(smech[1], props, dim)
                    Wi = S @ Wr
                    Xd, f = get_Xy(pmech[1], props, p_uf)
                    X = Xd @ Wi
                    NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wi
                if Type==3:
                    Wd = get_subspace(pmech[1], props, dim)
                    Wi = normalize(S.transpose() @ Wd)
                    Xd, f = get_Xy(pmech[1], props, p_uf)
                    X = Xd @ S @ Wi
                    NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
                if Type==4:
                    Wr = get_subspace(smech[1], props, dim)
                    Xr, f = get_Xy(smech[1], props, s_uf)
                    X = Xr @ S.transpose() @ S @ Wr
                    NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ S @ Wr
                try:
                    response_surface.train(X,f)
                    v,dv = response_surface.predict(NX)
                    minA, maxA = min(minA, np.mean(v)), max(maxA, np.mean(v))
                    minB, maxB = min(minB, np.std(v))*0.8, max(maxB, np.std(v)*1.2)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                   print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                   data_list.append([None,None])
                   continue
                dist_list.append([X[:,0],f])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[i][p].plot(T_revert, data_list[:,0],'C1')
            BX[i][p].plot(T_revert, data_list[:,1],'C1')
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],'C1.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])

    if Flag==True:
        Type = 2
        # = = = = = = = = = = = =
        # pmech -> intermediate 2 if necessary
        Legend.append("transition - "+TypeName[Type])
        cprint("\nPrinting transition - "+TypeName[Type],'g')
        data_dict, dist_dict = {}, {}
        for idx_phi,phi in enumerate(phi_arr):
            data_dict[phi], dist_dict[phi] = {}, {}

            for idx_P,P in enumerate(P_arr):
                data_list, dist_list = [], []

                for idx_T,T in enumerate(T_arr):
                    progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                    props['phi'],props['P'],props['T'] = phi,P,T
                    
                    if Type==1:
                        Wd = get_subspace(pmech[1], props, dim)
                        Wi = normalize(S.transpose() @ Wd)
                        Xr, f = get_Xy(smech[1], props, s_uf)
                        X = Xr @ Wi
                        NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
                    if Type==2:
                        Wr = get_subspace(smech[1], props, dim)
                        Wi = S @ Wr
                        Xd, f = get_Xy(pmech[1], props, p_uf)
                        X = Xd @ Wi
                        NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wi
                    if Type==3:
                        Wd = get_subspace(pmech[1], props, dim)
                        Wi = normalize(S.transpose() @ Wd)
                        Xd, f = get_Xy(pmech[1], props, p_uf)
                        X = Xd @ S @ Wi
                        NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
                    if Type==4:
                        Wr = get_subspace(smech[1], props, dim)
                        Xr, f = get_Xy(smech[1], props, s_uf)
                        X = Xr @ S.transpose() @ S @ Wr
                        NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ S @ Wr
                    try:
                        response_surface.train(X,f)
                        v,dv = response_surface.predict(NX)
                        minA, maxA = min(minA, np.mean(v)), max(maxA, np.mean(v))
                        minB, maxB = min(minB, np.std(v))*0.8, max(maxB, np.std(v)*1.2)
                        data_list.append([np.mean(v),np.std(v)])
                    except:
                       print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                       data_list.append([None,None])
                       continue
                    dist_list.append([X[:,0],f])
                data_dict[phi][P] = data_list
                dist_dict[phi][P] = dist_list

        # prepare figures
        T_revert = np.array([1000./Ti for Ti in T_arr])
        for i,phi in enumerate(phi_arr):
            for p,P in enumerate(P_arr):
                data_list = np.array(data_dict[phi][P])
                dist_list = dist_dict[phi][P]
                AX[i][p].plot(T_revert, data_list[:,0],'C3')
                BX[i][p].plot(T_revert, data_list[:,1],'C3')
                for t,T in enumerate(T_arr):
                    height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                    b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                    width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                    l_width = (0-min(dist_list[t][0]))/6*0.03
                    left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                    bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                    plt.figure("Data Distributions")
                    ax = plt.axes([left,bottom,width,height],facecolor='none')
                    plt.plot(dist_list[t][0],dist_list[t][1],'C3.',ms=1)
                    for edge in ['top','right','bottom','left']:
                        ax.spines[edge].set_visible(False)
                    plt.xticks([])
                    plt.yticks([])
        Type = 5

    # = = = = = = = = = = = =
    # smech
    Legend.append(smech[1])
    m,name,mech = smech
    cprint("\nPrinting %s"%name,'g')
    
    # Load data
    data_dict, dist_dict = {}, {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi], dist_dict[phi] = {}, {}
        
        for idx_P,P in enumerate(P_arr):
            data_list, dist_list = [], []

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
                    minA, maxA = min(minA, np.mean(v)), max(maxA, np.mean(v))
                    minB, maxB = min(minB, np.std(v))*0.8, max(maxB, np.std(v)*1.2)
                    data_list.append([np.mean(v),np.std(v)])
                except:
                    print("\033[1;35mError occured. Maybe not enought samples.\033[0m")
                    data_list.append([None,None])
                    continue
                dist_list.append([X[:,0],f])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[i][p].plot(T_revert, data_list[:,0],'C2')
            BX[i][p].plot(T_revert, data_list[:,1],'C2')
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],'C2.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])
    
    set_sub_plots(AX, r'1000/T, $K^{-1}$', r'$\log(\tau)$', Legend, xlim=[minT,maxT],ylim=[minA,maxA])
    set_sub_plots(BX, r'1000/T, $K^{-1}$', r'$\sigma of \log(\tau)$',Legend, xlim=[minT,maxT],ylim=[minB,maxB])
    save_figure(figA, figs_dir+'prop_%s_%s_Dist_%d.png'%(pmech[1],smech[1],Type))
    save_figure(figB, figs_dir+'prop_%s_%s_Sigm_%d.png'%(pmech[1],smech[1],Type))
    plt.show()
