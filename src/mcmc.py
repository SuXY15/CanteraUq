import pandas as pd
from utils import *
from respSurface import *

xm=np.array([float(v) for v in
    """
    -1.38740435 -0.33781194  0.55995156 -1.39986732 -0.33656576  0.85876974
    -0.26712536 -1.34194205 -0.37606529 -0.02585443 -0.24723983 -0.6508789
    -0.3844979   0.64331919  0.52430331 -0.21468851  0.51802463  0.9500797
    0.05225793  0.83298441  3.75166555  2.83092485  0.936531   -3.75589822
    3.71409842  1.17407511 -0.36412074  1.34771054  0.34404334 -0.53447187
    5.05730251 -0.11235795  2.00536432""".split()])

def singel_MCMC(props, dim, order, N, show=False):
    if show:
        fig,ax = plt.subplots(figsize=c2i(12,9))

    response_surface = PolynomialApproximation(N=order)
    cprint("Showing condition phi=%.1f p=%.1fatm T=%.1fK"%
            (props['phi'],props['P'],props['T']), 'g')
    
    # the true mech
    m,name,mech = mechs[0]
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%\
                (name,UF,props['phi'],props['P'],props['T'],samplingSize)
    uncetainty_factors = load_uncertainty(mech[:-3]+'txt')
    VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
    idt, factor, tdata = load_training_set(file_name)
    f = np.transpose([np.log10(idt),])
    X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
    response_surface.train(X, f)

    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
    NP, dn = response_surface.predict(NX)

    if show:
        dist = pd.DataFrame(NP.flatten().tolist(), columns=['H2'])
        dist.plot.kde(ax=ax, color='k', lw=1)


    # load current mech
    m,name,mech = mechs[1]
    file_name = acts_dir+"%s_UF=%.1f_phi=%.1f_p=%.1fatm_T=%.1fK_s=%d.json"%\
                (name,UF,props['phi'],props['P'],props['T'],samplingSize)
    uncetainty_factors = load_uncertainty(mech[:-3]+'txt')
    VT = np.array(json_reader(file_name[:-5]+"_as.json")[0]['subspace'])
    idt, factor, tdata = load_training_set(file_name)
    f = np.transpose([np.log10(idt),])
    X = np.dot([3.*np.log(fact)/np.log(uncetainty_factors) for fact in factor], np.transpose(VT[:dim]))
    response_surface.train(X, f)

    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
    NP, dn = response_surface.predict(NX)

    if show:
        dist = pd.DataFrame(NP.flatten().tolist(), columns=['H2fake'])
        dist.plot.kde(ax=ax, color='r', lw=1)

    plt.show()

if __name__ == "__main__":
    dim, order,N = 2,2,50000
    
    props['phi'], props['P'], props['T'] = 1.0, 10.0, 1200.

    singel_MCMC(props, dim, order, N, True)