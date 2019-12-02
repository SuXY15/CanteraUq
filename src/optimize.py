import cvxpy as cp
import pandas as pd
from qcqp import *
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

def opt_condtion(props, dim, order, N, mu_tol, sigma_tol, show=False):
    if show:
        fig,ax = plt.subplots(figsize=c2i(12,9))
    response_surface = PolynomialApproximation(N=order)
    cprint("Optimizing condition phi=%.1f p=%.1fatm T=%.1fK"%
            (props['phi'],props['P'],props['T']), 'g')
    
    # load target data
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

    mu_obs = np.mean(np.log10(idt))
    sigma_obs = np.std(np.log10(idt))

    print("target:\n\tmu: %.5e, sigma: %.5e"%(mu_obs, sigma_obs))

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
    
    y0 = response_surface.poly_weights[0][0]
    g = response_surface.g
    H = response_surface.H

    mu_cur = y0 + H.trace()
    sigma_cur = g.T @ g + 2*(H@H).trace()
    print("current:\n\tmu: %.5e, sigma: %.5e"%(mu_cur, sigma_cur))

    # # SOLVE BY CVXPY and QCQP
    # x = cp.Variable(dim)
    # c1 = y0 + H.trace() - mu_obs
    # c2 = g.T @ g + 2*(H@H).trace() - sigma_obs**2
    # cons = [cp.quad_form(x,H) + g.T * x + c1 <= mu_tol,
    #         cp.quad_form(x,H) + g.T * x + c1 >=-mu_tol,
    #         cp.quad_form(x,H@H)*4 + 4*g.T@H * x + c2 <= sigma_tol,
    #         cp.quad_form(x,H@H)*4 + 4*g.T@H * x + c2 >=-sigma_tol,
    #         ]
    # prob = cp.Problem(cp.Minimize(cp.sum_squares(x)), cons)

    # # Create a QCQP handler.
    # qcqp = QCQP(prob)
    # qcqp.suggest(SDR)
    
    # # Attempt to improve the starting point given by the suggest method
    # f_cd, v_cd = qcqp.improve(COORD_DESCENT)
    # # print("Coordinate descent: objective %.3f, violation %.3f" % (f_cd, v_cd))
    # # print("x.value:", x.value.flatten())
    # # print("solve x:", x.value.reshape(1,dim) @ VT[:dim])

    # # for solved x
    # x0 = x.value
    # SOLVE BY CVXPY and QCQP
    x = cp.Variable(dim+2)
    c1 = y0 + H.trace() - mu_obs
    c2 = g.T @ g + 2*(H@H).trace() - sigma_obs**2
    cons = [cp.quad_form(x[:dim],H) + g.T * x[:dim] + c1 <= x[dim],
            cp.quad_form(x[:dim],H) + g.T * x[:dim] + c1 >=-x[dim],
            cp.quad_form(x[:dim],H@H)*4 + 4*g.T@H * x[:dim] + c2 <= x[dim+1],
            cp.quad_form(x[:dim],H@H)*4 + 4*g.T@H * x[:dim] + c2 >=-x[dim+1],
            ]
    prob = cp.Problem(cp.Minimize(cp.sum_squares(x[dim:])/100+cp.sum_squares(x[:dim]))*100, cons)

    # Create a QCQP handler.
    qcqp = QCQP(prob)
    qcqp.suggest(SDR)
    
    # Attempt to improve the starting point given by the suggest method
    f_cd, v_cd = qcqp.improve(COORD_DESCENT)
    x0 = x.value[:dim]


    mu_opt = y0 + H.trace() + x0.T @ H @ x0 + g.T @ x0
    sigma_opt = np.sqrt(g.T @ g + 2*(H@H).trace() + 4*g.T @ H @ x0 + 4*x0.T @ H @ H @ x0)
    print("optimize error:")
    print("\tmu\t%.5e, %.5e, %.5e"%(mu_opt, mu_obs, mu_opt - mu_obs))
    print("\tsigma\t%.5e, %.5e, %.5e"%(sigma_opt, sigma_obs, sigma_opt**2 - sigma_obs**2))

    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
    NP = np.array([y0+(g+2*H@x0).T@ NXi.T + g.T @ x0 + NXi @ H @ NXi.T + x0.T @ H @ x0 for NXi in NX])

    if show:
        dist = pd.DataFrame(NP.flatten().tolist(), columns=['modified'])
        dist.plot.kde(ax=ax, color='b', lw=1)
    
    # for target x
    x0 = np.dot(VT[:dim],xm)
    mu_tgt = y0 + H.trace() + x0.T @ H @ x0 + g.T @ x0
    sigma_tgt = np.sqrt(g.T @ g + 2*(H@H).trace() + 4*g.T @ H @ x0 + 4*x0.T @ H @ H @ x0)
    print("target error:")
    print("\tmu\t%.5e, %.5e, %.5e"%(mu_tgt, mu_obs, mu_tgt - mu_obs))
    print("\tsigma\t%.5e, %.5e, %.5e"%(sigma_tgt, sigma_obs, sigma_tgt - sigma_obs))

    NX = np.transpose(np.dot(VT[:dim], np.random.normal(0.0, 1.0, (len(uncetainty_factors),N))))
    NP = np.array([y0+(g+2*H@x0).T@ NXi.T + g.T @ x0 + NXi @ H @ NXi.T + x0.T @ H @ x0 for NXi in NX])
    if show:
        dist = pd.DataFrame(NP.flatten().tolist(), columns=['real $x_0$'])
        dist.plot.kde(ax=ax, color='m', lw=1)
    
    if show:
        plt.text(-4.5, 1.8, r'$\phi=%.1f$'%props['phi'], fontsize=15)
        plt.text(-4.5, 1.4, r'$P=%datm$'%props['P'], fontsize=15)
        plt.text(-4.5, 1.0, r'$T_0=%dK$'%props['T'], fontsize=15)
        plt.xlim([-6,-3.5])

        # plt.text(-4.3, 3.8, r'$\phi=%.1f$'%props['phi'], fontsize=15)
        # plt.text(-4.3, 2.4, r'$P=%datm$'%props['P'], fontsize=15)
        # plt.text(-4.3, 1.0, r'$T_0=%dK$'%props['T'], fontsize=15)
        # plt.xlim([-4.6,-4.1])
        
        plt.xlabel(r"$\log_{10}(\rm{IDT}[s])$")
        plt.ylabel(r"PDF")
        save_figure(fig, figs_dir+"optimize_pdf_phi=%.1f p=%.1fatm T=%.1fK.png"%
                         (props['phi'],props['P'],props['T']))
        plt.show()
    
    return VT[:dim], x.value[:dim].reshape(dim,1)

def opt_conditions(dim=2, order=2, N=50000):
    mu_tol = 1e-4
    sigma_tol = 2e-2

    A = None
    b = None
    idx = None
    x = None
    
    # conditions = [ [0.5,  1.0, 1200.],
    #                [0.5, 10.0, 1200.],
    #                [1.5, 20.0, 1200.],
    #                [1.5, 10.0, 1400.] ]
    # for props['phi'],props['P'],props['T'] in conditions:

    for props['phi'] in [0.5, 1.0, 1.5]:
        for props['P'] in [1.0, 10.0, 20.0]:
            for props['T'] in [1200., 1400., 1600.]:
                mu_tol = 1e-4
                sigma_tol = 2e-2
                Ai, bi = opt_condtion(props, dim, order, N, mu_tol, sigma_tol)
                if A is None:
                    idx = [i for i,v in enumerate(Ai[0]) if abs(v)>0.]
                    x = Ai[0]
                    Ai = np.array([ai[idx] for ai in Ai])
                    A = Ai
                    b = bi
                else:
                    Ai = np.array([ai[idx] for ai in Ai])
                    A = np.vstack((A, Ai))
                    b = np.vstack((b, bi))
                cprint(bi, 'b')

    x_opt = np.linalg.lstsq(A, b, rcond=-1)[0]
    print(x_opt)

    x[idx] = x_opt.flatten()
    print(x)
    return x

if __name__ == "__main__":
    dim, order,N = 2,2,50000

    # props['phi'], props['P'], props['T'] = 1.0, 10.0, 1200.
    # opt_condtion(props, dim, order, N, 0, 0, True)

    opt_conditions(dim, order, N)