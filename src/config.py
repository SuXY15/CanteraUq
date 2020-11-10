#
# For mechanism configuration, should be run before DRG
#
from utils import *

if len(sys.argv)<2:
    cprint("Error! Config need a mechanism name")
    sys.exit(0)

mech = sys.argv[1]

if __name__ == "__main__":
    if mech == "H2":
        UF = 0.
        pdiff = 5e-3
        maxt = 1e-2          # max time for json writer
        samplingSize = 100  # sampling size

        mech_dir = "mech/H2/"
        comp_dir = "data/H2_comp/" # comparation data
        sens_dir = "data/H2_sens/" # sensitivity data
        acts_dir = "data/H2_acts/" # active subspace data
        resp_dir = "data/H2_resp/" # response surface data
        figs_dir = "data/H2_figs/" # figures
        # mechanisms
        mech_arr = ["H2", "H2fake"]
        
        # conditions
        phi_arr=[0.5, 1.0, 1.5]
        P_arr = [1.0, 10.0, 20.0]
        T_arr = [1200., 1300., 1400., 1500., 1600.]
        props = {'T': 1200.,                # Temperature
                 'P': 10.,                  # Pressure
                 'phi': 1.0,                # equivalence ratio
                 'air':'O2:1.0,N2:3.76',    # air components
                 'fuel':'H2:1.0',      # fuel components
                 'type':'UV',               # sim type
                 'tot': .1,                # total time
                 'UF': UF,
             }

    # = = = = = = = = = =
    # DMEzhao
    elif mech == "DME":
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
                 'UF': UF,
             }
             
    # = = = = = = = = = =
    # DME2000
    elif mech == "DME2000":
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
                    "DMEsk50",
                    "DMEsk43",
                    "DMEsk29"]
        # conditions
        phi_arr=[0.5, 1.0, 1.5]
        P_arr = [1.0, 5.0, 20.0]
        T_arr = [600., 650., 700., 800., 1000., 1200.]
        
        props = {'T': 1000.,                # Temperature
                 'P': 10.,                  # Pressure
                 'phi': 1.0,                # equivalence ratio
                 'air':'O2:1.0,N2:3.76',    # air components
                 'fuel':'CH3OCH3:1.0',      # fuel components
                 'type':'UV',               # sim type
                 'tot': 10.,                # total time
                 'UF': UF,
             }
    # = = = = = = = = = =
    # Methane without N2 branch
    elif mech == "CH4":
        UF = 2.             # Uncertainty factor
        pdiff = 5e-3        # difference interval for sensitivity
        maxt = 5e-2         # max time for json writer
        samplingSize = 400  # sampling size

        mech_dir = "mech/CH4/"
        comp_dir = "data/CH4_comp/" # comparation data
        sens_dir = "data/CH4_sens/" # sensitivity data
        acts_dir = "data/CH4_acts/" # active subspace data
        resp_dir = "data/CH4_resp/" # response surface data
        figs_dir = "data/CH4_figs/" # figures
        
        # mechanisms
        mech_arr = ["CH4_35",
                    "CH4_33",
                    "CH4_31",
                    "CH4_30",
                    "CH4_28",
                    "CH4_26",
                    "CH4_25",
                    "CH4_24",
                    "CH4_22",
                    "CH4_19",
                    "CH4_18",
                    "CH4_16",
                    "CH4_15"]

        # conditions
        phi_arr=[0.5, 1.0, 1.5]
        P_arr = [1.0, 10.0, 20.0]
        T_arr = [1000., 1100., 1200., 1300., 1400.]
        
        props = {'T': 1200.,                # Temperature
                 'P': 10.,                  # Pressure
                 'phi': 1.0,                # equivalence ratio
                 'air':'O2:1.0,N2:3.76',    # air components
                 'fuel':'CH4:1.0',      # fuel components
                 'type':'UV',               # sim type
                 'tot': 10.,                # total time
                 'UF': UF,
             }
    elif mech == "JP10":
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
    checkpath(comp_dir); checkpath(sens_dir)
    checkpath(acts_dir); checkpath(resp_dir)
    checkpath(figs_dir); checkpath(mech_dir)
    config_writer('data/'+mech+".config", config_dict)
    conditions_writer(mech_dir, props, get_conditions(phi_arr,P_arr,T_arr))
    print(json.dumps(config_reader('data/'+mech+".config"), indent=4))
