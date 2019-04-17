from utils import *

choice = sys.argv[1]

if __name__ == "__main__":
    if choice == "DME":
        UF = 5.             # Uncertainty factor
        pdiff = 5e-3        # difference interval for sensitivity
        maxt = 5e-2         # max time for json writer
        samplingSize = 400  # sampling size

        mech_dir = "mech/DME/"
        comp_dir = "data/DME_comp/" # comparation data
        sens_dir = "data/DME_sens/" # sensitivity data
        acts_dir = "data/DME_acts/" # active subspace data
        resp_dir = "data/DME_resp/" # response surface data
        figs_dir = "data/DME_figs/" # figures
        
        # mechanisms
        mech_arr = ["DMEzhao",
                    "DMEsk42",
                    "DMEsk40",
                    "DMEsk30"]
        # conditions
        phi_arr=[0.5, 1.0, 1.5]
        P_arr = [1.0, 10.0, 20.0]
        T_arr = [650., 700., 800., 1000., 1200.]
        
        props = {'T': 1000.,                # Temperature
                 'P': 10.,                  # Pressure
                 'phi': 1.0,                # equivalence ratio
                 'air':'O2:1.0,N2:3.76',    # air components
                 'fuel':'CH3OCH3:1.0',      # fuel components
                 'type':'UV',               # sim type
                 'tot': 10.,                # total time
             }
    elif choice == "DME2000":
        UF = 5.             # Uncertainty factor
        pdiff = 1e-2        # difference interval for sensitivity
        maxt = 5e-2         # max time for json writer
        samplingSize = 512  # sampling size

        mech_dir = "mech/DME2000/"
        comp_dir = "data/DME2000_comp/" # comparation data
        sens_dir = "data/DME2000_sens/" # sensitivity data
        acts_dir = "data/DME2000_acts/" # active subspace data
        resp_dir = "data/DME2000_resp/" # response surface data
        figs_dir = "data/DME2000_figs/" # figures
        
        # mechanisms
        mech_arr = ["DME2000",
                    "DMEsk59",
                    "DMEsk53",
                    "DMEsk47"]
        # conditions
        phi_arr=[0.5, 1.0, 1.5]
        P_arr = [1.0, 10.0, 20.0]
        T_arr = [650., 700., 800., 1000., 1200.]
        
        props = {'T': 1000.,                # Temperature
                 'P': 10.,                  # Pressure
                 'phi': 1.0,                # equivalence ratio
                 'air':'O2:1.0,N2:3.76',    # air components
                 'fuel':'CH3OCH3:1.0',      # fuel components
                 'type':'UV',               # sim type
                 'tot': 10.,                # total time
             }
    elif choice == "JP10":
        raise Exception("JP10 not defined.")
    else:
        raise Exception("No choices matched.")
        
    config_dict = { "UF": UF,
                    "pdiff": pdiff,
                    "maxt": maxt,
                    "samplingSize": samplingSize,
                    "comp_dir": comp_dir,
                    "sens_dir": sens_dir,
                    "acts_dir": acts_dir,
                    "resp_dir": resp_dir,
                    "figs_dir": figs_dir,
                    "mech_dir": mech_dir,
                    "mech_arr": mech_arr,
                    "phi_arr": phi_arr,
                    "P_arr": P_arr,
                    "T_arr": T_arr,
                    "props": props
                }

    config_writer(choice+".config",config_dict)
    conditions_writer(mech_dir,props,get_conditions(phi_arr,P_arr,T_arr),note='')
    print(json.dumps(config_reader(choice+".config"),indent=4))
