from utils import *
sys.path.append("dep/active_subspaces/utils")
from response_surfaces import *


if __name__ == "__main__":
    gas0 = load_mech(mechs[0][2])
    eqs0 = [r.equation for r in gas0.reactions()]
    main_reaction = set()
    
    for m,name,mech in mechs:
        print("Mech: %s"%mech)

        # loading
        gas = load_mech(mech)
        filename = sens_dir+name+".json"
        for data in json_reader(filename):
            for p in  [i for i,di in enumerate(np.isinf(data['tdata'])) if di]:
                props = data['props']
                idt = props['idt']
                Tcurv = getTcurv(gas, props)
                
                gas.set_multiplier(props['factor'][p]*(1+props['pdiff']), p)
                try:
                    pidt = get_ign(gas, props)
                    pTcurv = getTcurv(gas, props)
                    sens = (np.log(pidt) - np.log(idt))/np.log(1+props['pdiff'])
                except:
                    gas.set_multiplier(props['factor'][p], p)
                    gas.set_multiplier(props['factor'][p]*(1+props['pdiff']/2), p)
                    pidt = get_ign(gas, props)
                    sens = (np.log(pidt) - np.log(idt))/np.log(1+props['pdiff']/2)

                gas.set_multiplier(props['factor'][p], p)
                
                cprint("%s: phi %.1f P %.1f T %.1f"%(mech, props['phi'], props['P'], props['T']))
                cprint("oris %.4e idt %.6e, pidt %.6e, sens %.4e"%(data['tdata'][p], idt, pidt, sens))
                # plt.plot(Tcurv['t'], Tcurv['T'], 'g')
                # plt.plot(pTcurv['t'], pTcurv['T'],'r')

                # plt.show()

