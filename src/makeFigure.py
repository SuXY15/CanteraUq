from utils import *

legend_format = {'handletextpad':0.1, 'labelspacing':0.4, 'columnspacing':1.0}
matplotlib.rc('legend', **legend_format )

def get_sub_plots3(num):
    fig, axs = plt.subplots(len(P_arr), 1, figsize=c2i(12,15), num=num)
    fig.canvas.set_window_title(num)
    return fig, axs

def set_sub_plots2(axs, xlabel, ylabel, legend, xlim=None, ylim=None, xscale=None, \
                    yscale=None, loc="upper center", anchor=[0.5,0.9], ncol=2):
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            if(i!=2):     axs[i][p].get_xaxis().set_visible(False)
            if(p!=2):     axs[i][p].get_yaxis().set_visible(False)
            if(i==1 and p==2):      axs[i][p].set_ylabel(ylabel)
            if ylim:   axs[i][p].set_ylim(ylim)
            if yscale: axs[i][p].set_yscale(yscale)

def set_sub_plots3(axs, xlabel, ylabel, legend, xlim=None, ylim=None, xscale=None, \
                    yscale=None, loc="upper center", anchor=[0.5,0.9], ncol=2):
    """ Set subplots
    """
    i,phi = 0, phi_arr[1]
    for p,P in enumerate(P_arr):
        if(p!=2):     axs[p].get_xaxis().set_visible(False)
        if p==2:      axs[p].set_xlabel(xlabel)
        if p==1:      axs[p].set_ylabel(ylabel)
        if p==1:      axs[p].legend(legend, loc=loc, bbox_to_anchor=anchor, frameon=False, ncol=ncol)
        if xlim:   axs[p].set_xlim(xlim)
        if ylim:   axs[p].set_ylim(ylim)
        if xscale: axs[p].set_xscale(xscale)
        if yscale: axs[p].set_yscale(yscale)
        axs[p].set_title(r"$\phi=%.1f, P=%.0fatm$"%(phi,P), pad=-15)

# Not used
# Figure showing active subspace
def showAS():
    font={'size':15}
    matplotlib.rc('font', **font)
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x = np.arange(-1, 1, 0.1)
    y = np.arange(-1, 1, 0.1)

    X = np.zeros([len(x), len(y)])
    Y = np.zeros([len(x), len(y)])
    Z = np.zeros([len(x), len(y)])

    for i,xi in enumerate(x):
        for j,yj in enumerate(y):
            X[i,j] = xi;
            Y[i,j] = yj;
    Z = np.e**(X+Y)

    ax.plot_surface(X,Y,Z,color='k',alpha=0.5)
    ax.plot(x,y,np.e**(x+y), 'r-')
    ax.plot(x,y,x*0.0, 'r--')
    ax.plot([x[0],x[0]],[y[0],y[0]],[0, np.e**(x[0]+y[0])],'r--')
    ax.plot([x[-1],x[-1]],[y[-1],y[-1]],[0, np.e**(x[-1]+y[-1])],'r--')
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")
    ax.set_zlabel("$f(x_1,x_2)$")
    save_figure(fig, "figures/Fig0_showAS.png")
    plt.show()

def DRG_curv(idx=[33,34,73]):
    font={'size':15}
    matplotlib.rc('font', **font)

    with open(mech_dir+'curv.txt') as f:
        data = []
        for line in f.readlines():
            lin = line.split(' ')
            lin[-1] = lin[-1].split('\n')[0]
            data.append([float(l) for l in lin if l!='' and l[0]>='0' and l[0]<='9' ])
    data = np.array(data)
    data[:,2] /= 100

    fig = plt.figure(figsize=c2i(12,9))

    ax1 = fig.add_subplot(111)
    ax1.plot(data[:,0], data[:,1], 'k-')
    ax1.plot(data[idx,0],data[idx,1],'ks',ms=4)
    ax1.set_xlabel(r'$\varepsilon_{DRG}$')
    ax1.set_ylabel(r'$N_{S}$')

    ax2 = ax1.twinx()
    ax2.plot(data[:,0], data[:,2]*100, 'r--')
    ax2.plot(data[idx,0],data[idx,2]*100,'ro',ms=4, )#markerfacecolor='none')
    ax2.set_ylabel(r'$err_{\max}$($\%$)')
    ax2.set_yscale(r'log')

    ax1.set_ylim(data[-1,1]-2,data[0,1]+2)
    ax2.set_ylim(1e-2, 1e6)
    def y1v2y2(y1v):
        return 10**((y1v-data[-1,1]+2)/(data[0,1]-data[-1,1]+4)*8-4)

    # for idxi in idx:
    #     ax2.plot([data[idxi,0], data[idxi,0]], [data[idxi,2], y1v2y2(data[idxi,1])], 'k--')

    # ax1.text(data[32,0]-0.11, data[32,1]-2, 'DME42', fontsize=15)
    # ax1.text(data[33,0]+0.02, data[33,1]-3, 'DME40', fontsize=15)
    # ax1.text(data[72,0]+0.02, data[72,1]-1, 'DME30', fontsize=15)

    fig.subplots_adjust(left=0.14,bottom=0.15,top=0.90,right=0.84,hspace=0.,wspace=0.)
    save_figure(fig, path='./figures/Fig2_DRG.png')
    plt.show()

mech_names = ['DME55', 'DME42', 'DME40', 'DME30']
def Tcurv_mech(props, fig, ax1setting, ax2setting, figname='(a)'):
    ax1t, ax1T, ax1pos = ax1setting
    ax1 = fig.add_axes(ax1pos)
    if ax2setting is not None:
        ax2t, ax2T, ax2pos = ax2setting
        ax2 = fig.add_axes(ax2pos)
    legend = []
    for m,name,mech in mechs:
        gas = load_mech(mech)
        data = getTcurv(gas, props)
        pos,idt = diffMax(data['t'], data['T'])
        ax1.plot(data['t'], data['T'], color_arr[m]+line_arr[m])
        if ax2setting is not None:
            ax2.plot(data['t'], data['T'], color_arr[m]+line_arr[m])
    for m,name,mech in mechs:
        gas = load_mech(mech)
        data = getTcurv(gas, props)
        pos,idt = diffMax(data['t'], data['T'])
        ax1.plot(idt, data['T'][pos], color_arr[m]+symbol_arr[m], fillstyle='none')
        legend.append(mech_names[m])
    ax1.set_xlim(ax1t)
    ax1.set_ylim(ax1T)
    if ax2setting is not None:
        ax2.set_xlim(ax2t); # ax2.set_xticks(ax2t)
        ax2.set_ylim(ax2T); # ax2.set_yticks(ax2T)
        ax1.plot(ax2t,[ax2T[0],ax2T[0]],'k-',lw=0.5); ax1.plot([ax2t[0],ax2t[0]],ax2T,'k-',lw=0.5)
        ax1.plot(ax2t,[ax2T[1],ax2T[1]],'k-',lw=0.5); ax1.plot([ax2t[1],ax2t[1]],ax2T,'k-',lw=0.5)
    if ax2setting == None:
        ax1.legend(legend, frameon=False, loc='best')
    ax1.text(ax1t[1]*0.72, ax1T[0]+(ax1T[1]-ax1T[0])*0.28, r'$\phi=%d$'%(props['phi']), fontsize=15)
    ax1.text(ax1t[1]*0.72, ax1T[0]+(ax1T[1]-ax1T[0])*0.16, r'$P=%datm$'%(props['P']), fontsize=15)
    ax1.text(ax1t[1]*0.72, ax1T[0]+(ax1T[1]-ax1T[0])*0.04, r'$T_u=%dK$'%(props['T']), fontsize=15)
    ax1.text(ax1t[0]+(ax1t[1]-ax1t[0])*0.01, ax1T[0]+(ax1T[1]-ax1T[0])*0.90, figname, fontsize=15)
    ax1.set_yticks([1000, 2000, 3000])
    ax1.set_xlabel(r'$t (s)$')
    ax1.set_ylabel(r'$T (K)$')

def Pathway():
    font={'size':15}
    matplotlib.rc('font', **font)

    # Low temperature
    fig = plt.figure("Fig2a", figsize=c2i(12,6))
    # props['phi'],props['P'],props['T'] = 1.0, 1.0, 600.0
    # ax1t=[0.00,0.15]; ax1T=[500,3000]
    # ax2t=[0.09,0.10]; ax2T=[600,1000]
    # ax1pos=[0.10, 0.10, 0.80, 0.80]
    # ax2pos=[0.25, 0.25, 0.20, 0.50]
    props['phi'],props['P'],props['T'] = 1.0, 10.0, 650.0
    ax1t=[0.00,0.020]; ax1T=[600,3000]
    ax2t=[0.0117,0.0131]; ax2T=[650,1150]
    ax1pos=[0.10, 0.10, 0.80, 0.80]
    ax2pos=[0.25, 0.25, 0.20, 0.50]

    Tcurv_mech(props, fig, [ax1t, ax1T, ax1pos], [ax2t, ax2T, ax2pos], figname='(b)')
    #fig.subplots_adjust(left=0.07,bottom=0.20,top=0.95,right=0.99,hspace=0.,wspace=0.)
    save_figure(fig, 'figures/Fig3b_pathway.png')

    # High temperature
    fig = plt.figure("Fig2b", figsize=c2i(12,6))
    props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200.0
    ax1t=[0.000, 0.002]; ax1T=[1000, 3000]
    ax1pos=[0.10, 0.10, 0.80, 0.80]
    Tcurv_mech(props, fig, [ax1t, ax1T, ax1pos], None, figname='(a)')
    #fig.subplots_adjust(left=0.07,bottom=0.20,top=0.95,right=0.99,hspace=0.,wspace=0.)
    save_figure(fig, 'figures/Fig3a_pathway.png')

    plt.show()

def compare():
    # = = = = = = = = = =
    # showing results
    IDT_DATA = {}
    maxA ,minA = -1e10, 1e10
    figA, AX = get_sub_plots3(num = "IDT comparasion")

    for m,name,mech in mechs:
        # loading data
        props_arr = [d['props'] for d in json_reader(comp_dir+name+".json")]
        rela_err = []
        i,phi = 0,phi_arr[1]
        if m==0: IDT_DATA[phi] = {}
        for p,P in enumerate(P_arr):
            # load data in array T and sorted by 1000/T
            props_parr = [props for props in props_arr if props['P']==P and props['phi']==phi]
            props_parr = sorted(props_parr, key=lambda props:1000./props['T'])
            Temp_arr = np.array([props['T'] for props in props_parr if props['T'] in T_arr])
            idt_arr = np.log10([props['idt'] for props in props_parr if props['T'] in T_arr])

            # plot IDT
            AX[p].plot(1000./Temp_arr,idt_arr,color_arr[m]+line_arr[m]+symbol_arr[m])
            maxA,minA = max(max(idt_arr)+0.9,maxA),min(min(idt_arr)-0.1,minA)

    # = = = = = = = = = =
    # figure setting
    figA.subplots_adjust(left=0.15,bottom=0.09,top=0.98,right=0.95,hspace=0.,wspace=0.)
    set_sub_plots3(AX, xlabel=r'$1000/T (K^{-1}$)', ylabel=r'$\log_{10}(\rm{IDT}[s])$',
                    legend=mech_names,ylim=[minA,maxA])
    save_figure(figA, path='figures/Fig3_compare_IDT.png')
    plt.show()

def subspace():
    font={'size':15}
    matplotlib.rc('font', **font)

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
    tdata_arr = normalize(tdata_arr[:10])
    C = np.transpose(tdata_arr) @ np.array(tdata_arr)

    # SVD decomposing
    U,S,VT = np.linalg.svd(C)
    # S,VT = np.linalg.eig(C) # resutls are the same
    
    # dimensional projection
    w1x = [np.dot(VT[0], 3.*np.log(props['factor'])/np.log(uncetainty_factors)) for props in props_arr]
    w2x = [np.dot(VT[1], 3.*np.log(props['factor'])/np.log(uncetainty_factors)) for props in props_arr]

    # figure(1): eigenvalue
    graph_rank = np.arange(1,21)
    tick = np.arange(0,21,2)
    fig1 = plt.figure("Eigenvalue",figsize=c2i(12,6))
    plt.plot(graph_rank,S[graph_rank-1], marker='o', markerfacecolor='none', color='k')
    plt.xticks(tick)
    plt.yscale(r'log')
    plt.xlabel(r'index')
    plt.ylabel(r'eigenvalue')
    fig1.subplots_adjust(left=0.18,bottom=0.25,top=0.95,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig1, path='figures/Fig4a_eigenvalue.png')

    # figure(2): one dimensional projection
    fig2 = plt.figure(r"$\mathbf{w}_1$ projection",figsize=c2i(12,6))
    plt.scatter(w1x, np.log10(idt), marker='o', color='', edgecolors='k')
    plt.xlabel(r'$w_1^T {x}$')
    plt.ylabel(r'$\log_{10}(\rm{IDT}[s])$')
    fig2.subplots_adjust(left=0.18,bottom=0.25,top=0.95,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig2, path='figures/Fig4c_w1x.png')

    # figure(3): w1 components
    x,y = np.transpose([[i,v] for i,v in enumerate(VT[0]) if abs(v)>0.05])
    eqs = [r.equation for r in gas.reactions()]
    for i,eq in enumerate(eqs):
        print(i, VT[0][i], eq)

    eqs = [r.equation for r in gas.reactions()]
    sVT = sorted([(i,v) for i,v in enumerate(VT[0])], key=lambda vi:-abs(vi[1]))

    fig3 = plt.figure("w1 components",figsize=c2i(12,6))
    markerline, stemlines, baseline = plt.stem(x,y, markerfmt='ko', linefmt='k-.',basefmt='gray')
    plt.plot([1,gas.n_reactions],[0,0],'-',color='gray')
    plt.setp(markerline, color='k', markerfacecolor='none', linewidth=2)
    plt.xlabel(r'reaction index')
    plt.ylabel(r'$w_1$ components')
    # plt.ylim([-0.65,0.5])
    plt.yticks([-0.5,0.0,0.5])
    plt.xlim([1, gas.n_reactions])
    fig3.subplots_adjust(left=0.18,bottom=0.22,top=0.95,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig3, path='figures/Fig4b_w1.png')
    plt.show()

def pdfplot():
    font={'size':15}
    matplotlib.rc('font', **font)
    
    import pandas as pd
    dim, order, N = 3,2,50000
    ind = 100
    response_surface = PolynomialApproximation(N=order)

    props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200.0

    # preparing
    p_pos, s_pos = 1, 2
    pmech = mechs[p_pos]
    smech = mechs[s_pos]
    
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

    fig, ax = plt.subplots(figsize=c2i(12,6))

    # # loading pmech
    # Wd = get_subspace(pmech[1], props, dim)
    # Xd, f = get_Xy(pmech[1], props, p_uf)
    # X = Xd @ Wd
    # NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wd
    # response_surface.train(X, f)
    # v,dv = response_surface.predict(NX)

    # print("Sigma p:", np.std(v))
    # dist = pd.DataFrame(v-np.mean(v))
    # dist.plot.kde(ax=ax, legend=False, style='k^-', lw=1, ind=ind, fillstyle='none')
    # #dist.plot.hist(density=True, ax=ax, bins=32, color='k', histtype='step', legend=False)

    # # transition
    # Wd = get_subspace(pmech[1], props, dim)
    # Wi = S.transpose() @ Wd
    # Xd, f = get_Xy(pmech[1], props, p_uf)
    # X = Xd @ S @ Wi
    # NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
    # response_surface.train(X, f)
    # v,dv = response_surface.predict(NX)

    # print("Sigma t:", np.std(v))
    # dist = pd.DataFrame(v-np.mean(v))
    # dist.plot.kde(ax=ax, legend=False, style='ko--', lw=1, ind=ind, fillstyle='none')
    # #dist.plot.hist(density=True, ax=ax, bins=32, color='k', style='--', histtype='step', legend=False)

    # # loading smech
    # Ws = get_subspace(smech[1], props, dim)
    # Xs, f = get_Xy(smech[1], props, s_uf)
    # X = Xs @ Ws
    # NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Ws
    # response_surface.train(X, f)
    # v,dv = response_surface.predict(NX)

    # print("Sigma s:", np.std(v))
    # dist = pd.DataFrame(v-np.mean(v))
    # dist.plot.kde(ax=ax, legend=False, style='r^-', lw=1, ind=ind, fillstyle='none')
    # #dist.plot.hist(density=True, ax=ax, bins=32, color='r', histtype='step', legend=False)


    # = = = = = = = = = = = = = = = = = = = =
    props['phi'],props['P'],props['T'] = 1.0, 10.0, 650.0
    # loading pmech
    Wd = get_subspace(pmech[1], props, dim)
    Xd, f = get_Xy(pmech[1], props, p_uf)
    X = Xd @ Wd
    NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wd
    response_surface.train(X, f)
    v,dv = response_surface.predict(NX)

    print("Sigma p:", np.std(v))
    dist = pd.DataFrame(v-np.mean(v))
    dist.plot.kde(ax=ax, legend=False, style='k^-', lw=1, ind=ind, fillstyle='none')
    #dist.plot.hist(density=True, ax=ax, bins=32, color='k', histtype='step', legend=False)

    # transition
    Wd = get_subspace(pmech[1], props, dim)
    Wi = S.transpose() @ Wd
    Xd, f = get_Xy(pmech[1], props, p_uf)
    X = Xd @ S @ Wi
    NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi
    response_surface.train(X, f)
    v,dv = response_surface.predict(NX)

    print("Sigma t:", np.std(v))
    dist = pd.DataFrame(v-np.mean(v))
    dist.plot.kde(ax=ax, legend=False, style='ko--', lw=1, ind=ind, fillstyle='none')
    #dist.plot.hist(density=True, ax=ax, bins=32, color='k', style='--', histtype='step', legend=False)

    # loading smech
    Ws = get_subspace(smech[1], props, dim)
    Xs, f = get_Xy(smech[1], props, s_uf)
    X = Xs @ Ws
    NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Ws
    response_surface.train(X, f)
    v,dv = response_surface.predict(NX)

    print("Sigma s:", np.std(v))
    dist = pd.DataFrame(v-np.mean(v))
    dist.plot.kde(ax=ax, legend=False, style='r^-', lw=1, ind=ind, fillstyle='none')
    #dist.plot.hist(density=True, ax=ax, bins=32, color='r', histtype='step', legend=False)
    
    plt.xlabel(r'$\log_{10}(\rm{IDT}[s])-\left<\log_{10}(\rm{IDT}[s])\right>$')
    plt.ylabel(r'PDF')
    plt.xlim([-0.5, 0.5])
    plt.ylim([-0.2, 3.5])
    # plt.legend([mech_names[p_pos], 'transition', mech_names[s_pos]], frameon=False)

    # plt.text(-0.48, 3.1, r'$(a)$', fontsize=15)
    # plt.text(-0.48, 2.6, r'$\phi=1$', fontsize=15)
    # plt.text(-0.48, 2.15, r'$P=1atm$', fontsize=15)
    # plt.text(-0.48, 1.7, r'$T_u=1200K$', fontsize=15)

    plt.text(-0.48, 3.1, r'$(b)$', fontsize=15)
    plt.text(-0.48, 2.6, r'$\phi=1$', fontsize=15)
    plt.text(-0.48, 2.15, r'$P=10atm$', fontsize=15)
    plt.text(-0.48, 1.7, r'$T_u=650K$', fontsize=15)

    fig.subplots_adjust(left=0.11,bottom=0.23,top=0.95,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig, path='figures/Fig8b_distribution.png')

    plt.show()

from respSurface import *
def trainPredict(dim=3,order=2,N=50000,method=''):
    font={'size':15}
    matplotlib.rc('font', **font)

    m,name,mech = mechs[0]
    props = config_dict['props']
    props['phi'], props['P'], props['T'] = 1.0, 10.0, 1200.
    response_surface = PolynomialApproximation(N=order)

    # prepare data
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%(name,
                UF,props['phi'],props['P'],props['T'],samplingSize)
    uncetainty_factors = load_uncertainty(mech, UF=UF)
    VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
    idt, factor, tdata = load_training_set(file_name)

    f = np.transpose([np.log10(idt),])
    X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
    trainSize = len(f)>>1
    x,v = X[:trainSize,:],X[trainSize:,:]
    xf,vf = f[:trainSize,:],f[trainSize:,:]
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

    # plt.xlim([-2.3, -0.6])
    fig1.subplots_adjust(left=0.18,bottom=0.22,top=0.95,right=0.95,hspace=0.,wspace=0.)
    fig2.subplots_adjust(left=0.18,bottom=0.22,top=0.95,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig1, "figures/M_Fig1_train&valid.png")
    save_figure(fig2, "figures/M_Fig1_prediction.png")
    plt.show()

def propagation():
    Legend = []
    dim, order, N = 3,2,50000
    font={'size':15}
    matplotlib.rc('font', **font)

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

    response_surface = PolynomialApproximation(N=order)

    # = = = = = = = = = = = =
    # pmech
    Legend.append(mech_names[p_pos])
    m,name,mech = pmech
    cprint("Printing %s"%name,'g')

    # Load data
    data_dict = {}
    maxB ,minB = -1e10, 1e10
    
    idx_phi, phi=1, phi_arr[1]
    data_dict[phi] = {}
    for idx_P,P in enumerate(P_arr):
        data_list = []
        for idx_T,T in enumerate(T_arr):
            props['phi'],props['P'],props['T'] = phi,P,T

            Wd = get_subspace(pmech[1], props, dim)
            Xd, f = get_Xy(pmech[1], props, p_uf)
            X = Xd @ Wd
            NX = np.transpose(np.random.normal(0.0, 1.0, (len(p_uf), N))) @ Wd

            response_surface.train(X, f)

            v,dv = response_surface.predict(NX)
            mean_v, std_v = response_surface.mean(Wd), response_surface.std(Wd)

            minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
            data_list.append([mean_v, std_v, np.mean(v), np.std(v)])
        data_dict[phi][P] = data_list
        
    # prepare figures
    figB, BX = get_sub_plots3(num="Data Sigmas")
    EX = []
    for p,P in enumerate(P_arr):
        EX.append(BX[p].twinx())
    
    minT, maxT = 1000./np.max(T_arr), 1000./np.min(T_arr)
    minT, maxT = maxT-(maxT-minT)*19/18, (maxT-minT)*19/18+minT
    T_revert = np.array([1000./Ti for Ti in T_arr])

    i, phi=1, phi_arr[1]
    for p,P in enumerate(P_arr):
        data_list = np.array(data_dict[phi][P])
        BX[p].plot(T_revert, data_list[:,1],color_arr[0])
    P_DATA_DICT = deepcopy(data_dict)

    # = = = = = = = = = = = =
    # pmech -> intermediate
    Legend.append("transition")
    cprint("\nPrinting transition",'g')
    data_dict = {}

    idx_phi, phi=1, phi_arr[1]    
    data_dict[phi] = {}
    for idx_P,P in enumerate(P_arr):
        data_list = []
        for idx_T,T in enumerate(T_arr):
            props['phi'],props['P'],props['T'] = phi,P,T

            Wd = get_subspace(pmech[1], props, dim)
            Wi = S.transpose() @ Wd
            Xd, f = get_Xy(pmech[1], props, p_uf)
            X = Xd @ S @ Wi
            NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wi

            response_surface.train(Xd @ Wd, f)

            v,dv = response_surface.predict(NX)
            mean_v, std_v = response_surface.mean(S @ Wi), response_surface.std(S @ Wi)

            minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
            data_list.append([mean_v, std_v, np.mean(v), np.std(v)])
        
        data_dict[phi][P] = data_list

    I_DATA_DICT = deepcopy(data_dict)

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    i, phi = 1, phi_arr[1]
    for p,P in enumerate(P_arr):
        data_list = np.array(data_dict[phi][P])
        BX[p].plot(T_revert, data_list[:,1],color_arr[0]+'--')

    # = = = = = = = = = = = =
    # smech
    Legend.append(mech_names[s_pos])
    m,name,mech = smech
    cprint("\nPrinting %s"%name,'g')

    # Load data
    data_dict = {}
    
    idx_phi, phi=1, phi_arr[1]
    data_dict[phi] = {}
    for idx_P,P in enumerate(P_arr):
        data_list = []
        for idx_T,T in enumerate(T_arr):
            props['phi'],props['P'],props['T'] = phi,P,T

            Wr = get_subspace(smech[1], props, dim)
            Xr, f = get_Xy(smech[1], props, s_uf)
            X = Xr @ Wr
            NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wr

            response_surface.train(X, f)

            v,dv = response_surface.predict(NX)
            mean_v, std_v = response_surface.mean(Wr), response_surface.std(Wr)

            minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
            data_list.append([mean_v, std_v, np.mean(v), np.std(v)])

        data_dict[phi][P] = data_list
    
    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    i,phi = 1, phi_arr[1]
    for p,P in enumerate(P_arr):
        data_list = np.array(data_dict[phi][P])
        BX[p].plot(T_revert, data_list[:,1],color_arr[1])

    S_DATA_DICT = deepcopy(data_dict)

    # set twinx 
    minE, maxE= 1e10, -1e10
    T_revert = np.array([1000./Ti for Ti in T_arr])
    i,phi=1, phi_arr[1]
    for p,P in enumerate(P_arr):
        p_data_list = np.array(P_DATA_DICT[phi][P])[:,1]
        i_data_list = np.array(I_DATA_DICT[phi][P])[:,1]
        s_data_list = np.array(S_DATA_DICT[phi][P])[:,1]
        c_data_list = np.abs(p_data_list-i_data_list)/(np.abs(p_data_list-i_data_list)+np.abs(i_data_list-s_data_list))
        #if(p==1): c_data_list[1] = np.nan
        EX[p].plot(T_revert, c_data_list, color_arr[1]+'^-')
        minE, maxE = min(minE, np.min(c_data_list))-0.02, max(maxE, np.max(c_data_list)*1.2)+0.02
    
    for p,P in enumerate(P_arr):
        EX[p].get_xaxis().set_visible(False)
        if p==1:   EX[p].set_ylabel('$r_t$')
        EX[p].set_xlim([minT, maxT])
        EX[p].set_ylim([minE, maxE])

    cprint("Setting plots", 'g')
    Legend = [r'$\sigma_{r,d}$', r'$\sigma_{r,t}$', r'$\sigma_{r,s}$']
    set_sub_plots3(BX, r'$1000/T, K^{-1}$', r'$\sigma_r$',Legend, xlim=[minT,maxT],ylim=[minB,maxB],ncol=3)
    cprint("Saving figures", 'g')
    figB.subplots_adjust(left=0.14, bottom=0.10, top=0.98, right=0.85, hspace=0., wspace=0.)
    save_figure(figB, 'figures/Fig6_%s_%s_Sigm.png'%(pmech[1],smech[1]))
    plt.show()

def propagation3():
    Legend = []
    dim, order, N = 3,2,50000
    font={'size':10}
    matplotlib.rc('font', **font)

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
                mean_v, std_v = response_surface.mean(Wd), response_surface.std(Wd)

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])

                v,dv = response_surface.predict(X)
                dist_list.append([list(f),list(v)])

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
            BX[p][i].plot(T_revert, data_list[:,1],color_arr[0])
            AX[p][i].plot(T_revert, data_list[:,0],color_arr[0])
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = float(max(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                half_axes_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                half_axes_width = float(np.mean(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                left = i*0.3 + 0.03 + (1000./T-minT)/(maxT-minT)*0.3-half_axes_width
                bottom = (2-p)*0.3 + 0.03 + (data_list[t,0]-minA)/(maxA-minA)*0.3-half_axes_height

                plt.figure(num="Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[0]+'.',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])
    P_DATA_DICT = deepcopy(data_dict)

    font={'size':15}
    matplotlib.rc('font', **font)
    # = = = = = = = = = = = =
    # pmech -> intermediate
    Legend.append("transition")
    cprint("\nPrinting transition",'g')
    figE = plt.figure(figsize=c2i(12,9),num="Elimination")

    EX = []
    for i,phi in enumerate(phi_arr):
        ex = []
        for p,P in enumerate(P_arr):
            ex.append(BX[i][p].twinx())
        EX.append(ex)

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
                mean_v, std_v = response_surface.mean(S @ Wi), response_surface.std(S @ Wi)

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])

                r1 = np.linalg.norm(Wi[:,0])
                r2 = std_v/P_DATA_DICT[phi][P][idx_T][1]
                r3 = np.std(v)/P_DATA_DICT[phi][P][idx_T][3]
                #plt.scatter(r1, r3, marker='v', color='', edgecolors='k')
                plt.scatter(r1, r2, marker='^', color='', edgecolors='k')
                
                v,dv = response_surface.predict(X)
                dist_list.append([list(f),list(v)])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list
            print("%.1f %2.0f %.5f %.5f %.5f %.5f %.5f"%(phi,P,data_list[0][1],data_list[1][1],
                data_list[2][1],data_list[3][1],data_list[4][1]))

    #plt.legend(["MC", "PCE"], frameon=False)
    plt.xlabel(r"$\|P^T \mathbf{w}_{d,1}\|$")
    plt.ylabel(r"$\sigma_{r,t}/\sigma_{r,d}$")

    I_DATA_DICT = deepcopy(data_dict)

    font={'size':10}
    matplotlib.rc('font', **font)
    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[p][i].plot(T_revert, data_list[:,0],color_arr[0]+'--')
            BX[p][i].plot(T_revert, data_list[:,1],color_arr[0]+'--')
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = float(max(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                half_axes_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                half_axes_width = float(np.mean(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                left = i*0.3 + 0.03 + (1000./T-minT)/(maxT-minT)*0.3-half_axes_width
                bottom = (2-p)*0.3 + 0.03 + (data_list[t,0]-minA)/(maxA-minA)*0.3-half_axes_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[1]+'^',ms=1,fillstyle="none")
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])

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
                # progress((idx_T+idx_P*len(T_arr)+idx_phi*len(T_arr)*len(P_arr))/(len(phi_arr)*len(P_arr)*len(T_arr)))
                props['phi'],props['P'],props['T'] = phi,P,T

                Wr = get_subspace(smech[1], props, dim)
                Xr, f = get_Xy(smech[1], props, s_uf)
                X = Xr @ Wr
                NX = np.transpose(np.random.normal(0.0, 1.0, (len(s_uf), N))) @ Wr

                response_surface.train(X, f)

                v,dv = response_surface.predict(NX)
                mean_v, std_v = response_surface.mean(Wr), response_surface.std(Wr)

                minA, maxA = min(minA, mean_v), max(maxA, mean_v)
                minB, maxB = min(minB, std_v)*0.8, max(maxB, std_v*1.2)
                data_list.append([mean_v, std_v, np.mean(v), np.std(v)])
                
                v,dv = response_surface.predict(X)
                dist_list.append([list(f),list(v)])
            data_dict[phi][P] = data_list
            dist_dict[phi][P] = dist_list
            print("%.1f %2.0f %.5f %.5f %.5f %.5f %.5f"%(phi,P,data_list[0][1],data_list[1][1],
                data_list[2][1],data_list[3][1],data_list[4][1]))

    S_DATA_DICT = deepcopy(data_dict)

    # prepare figures
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            data_list = np.array(data_dict[phi][P])
            dist_list = dist_dict[phi][P]
            AX[p][i].plot(T_revert, data_list[:,0],color_arr[1])
            BX[p][i].plot(T_revert, data_list[:,1],color_arr[1])
            for t,T in enumerate(T_arr):
                height = float((max(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                width = float(max(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                half_axes_height = float((np.mean(dist_list[t][1])-min(dist_list[t][1]))/(maxA-minA)*0.3)
                half_axes_width = float(np.mean(dist_list[t][0])-min(dist_list[t][0]))/1*0.03
                left = i*0.3 + 0.03 + (1000./T-minT)/(maxT-minT)*0.3-half_axes_width
                bottom = (2-p)*0.3 + 0.03 + (data_list[t,0]-minA)/(maxA-minA)*0.3-half_axes_height
                plt.figure("Data Distributions")
                ax = plt.axes([left,bottom,width,height],facecolor='none')
                plt.plot(dist_list[t][0],dist_list[t][1],color_arr[2]+'v',ms=1)
                for edge in ['top','right','bottom','left']:
                    ax.spines[edge].set_visible(False)
                plt.xticks([])
                plt.yticks([])

    minE, maxE= 1e10, -1e10
    T_revert = np.array([1000./Ti for Ti in T_arr])
    for i,phi in enumerate(phi_arr):
        for p,P in enumerate(P_arr):
            p_data_list = np.array(P_DATA_DICT[phi][P])[:,1]
            i_data_list = np.array(I_DATA_DICT[phi][P])[:,1]
            s_data_list = np.array(S_DATA_DICT[phi][P])[:,1]
            t_data_list = np.abs(p_data_list-i_data_list)
            c_data_list = (p_data_list-i_data_list)/(np.abs(p_data_list-i_data_list)+np.abs(i_data_list-s_data_list))
            #EX[i][p].plot(T_revert, t_data_list, color_arr[0]+'v-')
            EX[p][i].plot(T_revert, c_data_list, color_arr[1]+'^-')
            #minE, maxE = min(min(minE, np.min(t_data_list)*0.8),np.min(t_data_list)*1.2), max(maxE, np.max(t_data_list)*1.2)
            minE, maxE = min(minE, np.min(c_data_list)-0.1), max(maxE, np.max(c_data_list)+0.1)

    font = {'size':10}
    matplotlib.rc('font', **font)
    set_sub_plots2(EX, r'1000/T, $K^{-1}$', r'$r_t$', Legend, xlim=[minT, maxT], ylim=[minE, maxE])
    cprint("Setting plots", 'g')
    set_sub_plots(AX, r'1000/T, $K^{-1}$', r'$\log(IDT[s])$', Legend, xlim=[minT,maxT],ylim=[minA,maxA])
    set_sub_plots(BX, r'1000/T, $K^{-1}$', r'$\sigma_r$',Legend, xlim=[minT,maxT],ylim=[minB,maxB])
    cprint("Saving figures", 'g')
    figB.subplots_adjust(left=0.06, bottom=0.09, top=0.98, right=0.90, hspace=0., wspace=0.)
    save_figure(figA,'figures/Fig7_prop_%s_%s_Dist.png'%(pmech[1],smech[1]))
    save_figure(figB, 'figures/Fig7_%s_%s_Sigm.png'%(pmech[1],smech[1]))
    save_figure(figE, 'figures/Fig8_%s_%s_Elim.png'%(pmech[1],smech[1]))
    plt.show()

pdfplot()