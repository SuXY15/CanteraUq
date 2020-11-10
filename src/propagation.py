from utils import *
sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *

dim, order, N = 3,2,50000

Legend = []
TypeName = ['O','I', 'II', 'III', 'IV']

if __name__ == "__main__":
    p_pos, s_pos = 1, 2
    pmech = mechs[p_pos]
    smech = mechs[s_pos]

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
                # progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wd = get_subspace(pmech[1], props, dim)
                Xd, f = get_Xy(pmech[1], props, p_uf)
                X = Xd @ Wd
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wd

                response_surface.train(X, f)
                
                v,dv = response_surface.predict(NX)
                # mean_v, std_v = np.mean(v), np.std(v)
                
                # print()
                mean_v, std_v = response_surface.mean(Wd), response_surface.std(Wd)
                # print("%.4f %.4f"%((mean_v-np.mean(v))/mean_v, (std_v-np.std(v))/std_v))

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])
                
                dist_list.append([NX[:,0],v])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list
            print("%.1f %2.0f %.5f %.5f %.5f %.5f %.5f"%(phi,P,data_list[0][1],data_list[1][1],
                data_list[2][1],data_list[3][1],data_list[4][1]))

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
            BX[i][p].plot(T_revert, data_list[:,1],color_arr[0])
            AX[i][p].plot(T_revert, data_list[:,0],color_arr[0])
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure(num="Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[0]+'.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])
    P_DATA_DICT = deepcopy(data_dict)

    # = = = = = = = = = = = =
    # pmech -> intermediate
    Legend.append("transition")
    cprint("\nPrinting transition",'g')
    figE = plt.figure(figsize=c2i(12,9),num="Elimination")

    data_dict, dist_dict = {}, {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi], dist_dict[phi] = {}, {}

        for idx_P,P in enumerate(P_arr):
            data_list, dist_list = [], []

            for idx_T,T in enumerate(T_arr):
                # progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wd = get_subspace(pmech[1], props, dim)
                Wi = S.transpose() @ Wd
                Xd, f = get_Xy(pmech[1], props, p_uf)
                X = Xd @ S @ Wi
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi

                response_surface.train(Xd @ Wd, f)

                v,dv = response_surface.predict(NX)
                
                # mean_v, std_v = np.mean(v), np.std(v)
                # print()
                mean_v, std_v = response_surface.mean(S @ Wi), response_surface.std(S @ Wi)
                # print("%.4f %.4f"%((mean_v-np.mean(v))/mean_v, (std_v-np.std(v))/std_v))

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])

                r1 = np.linalg.norm(Wi[:,0])
                r2 = std_v/P_DATA_DICT[phi][P][idx_T][1]
                r3 = np.std(v)/P_DATA_DICT[phi][P][idx_T][3]
                plt.scatter(r1, r2, marker='o', color='', edgecolors='r')
                plt.scatter(r1, r3, marker='o', color='', edgecolors='k')
                # print("%.1f %2.f %4.f %.5f %.5f %.4e"%(phi, P, T, r1, r2, r2-r1))

                dist_list.append([NX[:,0],v])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list
            print("%.1f %2.0f %.5f %.5f %.5f %.5f %.5f"%(phi,P,data_list[0][1],data_list[1][1],
                data_list[2][1],data_list[3][1],data_list[4][1]))

    plt.xlabel(r"$\|P \mathbf{w}_{d,1}\|$")
    plt.ylabel(r"$\sigma_{r,t}/\sigma_{r,d}$")

    I_DATA_DICT = deepcopy(data_dict)

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[i][p].plot(T_revert, data_list[:,0],color_arr[0]+'--')
            BX[i][p].plot(T_revert, data_list[:,1],color_arr[0]+'--')
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[1]+'o',ms=1,fillstyle="none")
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])

    # = = = = = = = = = = = =
    # smech
    Legend.append(smech[1])
    m,name,mech = smech
    cprint("\nPrinting %s"%name,'g')
    figC = plt.figure(figsize=c2i(12,9),num="Coupling")
    # Load data
    data_dict, dist_dict = {}, {}
    for idx_phi,phi in enumerate(phi_arr):
        data_dict[phi], dist_dict[phi] = {}, {}
        
        for idx_P,P in enumerate(P_arr):
            data_list, dist_list = [], []

            for idx_T,T in enumerate(T_arr):
                # progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T
                
                Wr = get_subspace(smech[1], props, dim)
                Xr, f = get_Xy(smech[1], props, s_uf)
                X = Xr @ Wr
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wr

                response_surface.train(X, f)
                
                v,dv = response_surface.predict(NX)

                # print()
                mean_v, std_v = response_surface.mean(Wr), response_surface.std(Wr)
                # print("%.4f %.4f"%((mean_v-np.mean(v))/mean_v, (std_v-np.std(v))/std_v))

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])

                Wd = get_subspace(pmech[1], props, dim)
                Wi = S.transpose() @ Wd
                r1 = np.sqrt(np.dot(Wi[:,0],Wr[:,0]))
                r2 = std_v/I_DATA_DICT[phi][P][idx_T][1]
                
                plt.plot(r1, r2, color_arr[idx_T]+'.')
                # print("%.1f %2.f %4.f %.5f %.5f %.4e"%(phi, P, T, r1, r2, r1-r2))

                dist_list.append([NX[:,0],v])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list
            print("%.1f %2.0f %.5f %.5f %.5f %.5f %.5f"%(phi,P,data_list[0][1],data_list[1][1],
                data_list[2][1],data_list[3][1],data_list[4][1]))

    plt.xlabel(r"$\sqrt{P^T \mathbf{w}_{d,1}\cdot \mathbf{w}_{s,1}}$")
    plt.ylabel(r"$\sigma_{r,s}/\sigma_{r,t}$")
    S_DATA_DICT = deepcopy(data_dict)

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[i][p].plot(T_revert, data_list[:,0],color_arr[1])
            BX[i][p].plot(T_revert, data_list[:,1],color_arr[1])
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                b_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = (max(dist_list[t][0])-min(dist_list[t][0]))/6*0.03
                l_width = (0-min(dist_list[t][0]))/6*0.03
                left = p*0.3 + 0.05 + (1000./T-minT)/(maxT-minT)*0.3-l_width
                bottom = (2-i)*0.3 + 0.05 + (data_list[t,0]-minA)/(maxA-minA)*0.3-b_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[1]+'.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])
    
    print()
    cprint("Setting plots", 'g')
    set_sub_plots(AX, r'1000/T, $K^{-1}$', r'$\log(IDT[s])$', Legend, xlim=[minT,maxT],ylim=[minA,maxA])
    set_sub_plots(BX, r'1000/T, $K^{-1}$', r'$\sigma_r$',Legend, xlim=[minT,maxT],ylim=[minB,maxB])
    cprint("Saving figures", 'g')
    save_figure(figA, figs_dir+'prop_%s_%s_Dist.png'%(pmech[1],smech[1]))
    save_figure(figB, figs_dir+'prop_%s_%s_Sigm.eps'%(pmech[1],smech[1]))
    save_figure(figC, figs_dir+'prop_%s_%s_Coup.eps'%(pmech[1],smech[1]))
    save_figure(figE, figs_dir+'prop_%s_%s_Elim.eps'%(pmech[1],smech[1]))
    plt.show()
