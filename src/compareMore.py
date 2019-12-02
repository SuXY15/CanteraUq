from utils import *
from check import *
from respSurface import *

def load_idts(file_name):
    sens_data_arr = json_reader(file_name)
    props_arr = [sd['props'] for sd in sens_data_arr]
    idx = [i for i,props in enumerate(props_arr) if props['idt']>0]
    idt = np.array([props_arr[i]['idt'] for i in idx])
    return idt

if __name__=="__main__":
    # = = = = = = = = = =
    # calculating
    figA, AX = get_sub_plots(num="IDT Compare")
    T_revert = np.array([1000./Ti for Ti in T_arr])
    Legend = []
    dim, order,N = 2,2,10000
    response_surface = PolynomialApproximation(N=order)
    maxA ,minA = -1e10, 1e10
    for m,name,mech in mechs:
        cprint("Loading mech:"+name, 'g')
        for p,P in enumerate(P_arr):
            for i,phi in enumerate(phi_arr):
                data_list = []
                for idx_T,T in enumerate(T_arr):
                    props['phi'],props['P'],props['T'] = phi,P,T
                    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                                    UF,props['phi'],props['P'],props['T'],samplingSize)
                    
                    uncetainty_factors = load_uncertainty(mech[:-3]+'txt')
                    VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
                    idt, factor, tdata = load_training_set(file_name)
                    f = np.transpose([np.log10(idt),])
                    X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
                    response_surface.train(X, f)

                    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
                    if m==0:
                        NP, dn = response_surface.predict(NX)  
                    else:
                        # optimize under T=1200; P=10; phi=1.0; 
                        # x0p = [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         -8.04644108e-01, -2.90735436e-17, -9.79817126e-02,  3.83095262e-02,
                        #         4.73202206e-21, -5.63087895e-03,  1.24627399e-01,  1.26349662e+00,
                        #         7.33505196e-02,  0.00000000e+00,  0.00000000e+00, -3.61153041e-01,
                        #         0.00000000e+00,  1.03016731e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,
                        #         -2.90387697e-02,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00]
                        
                        # optimize under T=1200, 1600; P=1, 20; phi=0.5, 1.0, 1.5
                        # x0p = [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #        -3.37481084e+00,  0.00000000e+00,  7.49651408e+00,  5.70080475e+00,
                        #        -1.42109554e-14,  3.25445458e+00, -5.05085947e+00,  8.00022663e-01,
                        #        -1.31382238e+00,  0.00000000e+00,  0.00000000e+00, -9.13635052e+00,
                        #         0.00000000e+00, -3.27245578e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         4.79425339e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00]

                        # test for opt error
                        x0p = [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                                0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                                2.90551911e-01,  0.00000000e+00, -6.27201525e+00, -9.16347622e+00,
                                3.55269854e-14,  2.14862782e+01, -4.61089023e-01,  7.43212437e-01,
                               -1.03512344e+01,  0.00000000e+00,  0.00000000e+00, -5.02015359e-01,
                                0.00000000e+00,  1.05214166e+00,  0.00000000e+00,  0.00000000e+00,
                                0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                                1.03748235e+01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                                0.00000000e+00]
                        
                        # # optimize under T=1200, 1600; P=1, 10, 20; phi=0.5, 1.0, 1.5
                        # x0p = [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         -3.41366633e+0,  0.00000000e+00,  6.18997425e+00,  3.80778007e+00,
                        #         2.66452876e-14,  2.92049205e+00, -3.13791257e+00,  3.94087570e-01,
                        #         -7.90225605e-1,  0.00000000e+00,  0.00000000e+00, -6.91602198e+00,
                        #         0.00000000e+00, -2.68938314e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         2.75434593e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                        #         0.00000000e+00,]
                        x0 = VT[:dim] @ x0p
                        y0 = response_surface.poly_weights[0][0]
                        g = response_surface.g
                        H = response_surface.H
                        NP = np.array([y0+(g+2*H@x0).T@ NXi.T + g.T @ x0 + NXi @ H @ NXi.T + x0.T @ H @ x0 for NXi in NX])
                    NP = NP.flatten()
                    mu, sigma = np.mean(NP), np.std(NP)

                    # idts = np.log10(load_idts(file_name))
                    # mu, sigma = np.mean(idts), np.std(idts)

                    data_list.append([mu, sigma])
                data_list = np.array(data_list)
                idt = data_list[:,0]
                sig = data_list[:,1]
                AX[p,i].plot(T_revert, idt, color_arr[m]+line_arr[m]+symbol_arr[m])
                AX[p,i].fill_between(T_revert, idt-sig, idt+sig, facecolor=color_arr[m],
                                        alpha=0.4, interpolate=True)
                maxA,minA = max(max(idt+sig)+0.1,maxA),min(min(idt-sig)-0.1,minA)
        Legend.append(name)

    set_sub_plots(AX, xlabel=r'$1000/T$ (K$^{-1}$)', ylabel=r'$\log_{10}(\rm{IDT}[s])$',
                    legend=Legend,ylim=[minA,maxA])
    figA.subplots_adjust(left=0.07,bottom=0.09,top=0.98,right=0.95,hspace=0.,wspace=0.)
    save_figure(figA, path=figs_dir+'compareMore_IDT.png')
    plt.show()
