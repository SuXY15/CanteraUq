#
# For comparation of detailed mechanism and skeletal mechanism
#
from utils import *
from check import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# = = = = = = = = = =
# calculator
def calculator(name, gas, props_arr, file_name):
    data_dict = {'props':[],'tdata':[]}
    NUM = len(props_arr)
    for i,props in enumerate(props_arr):
        print("|{0:07.3f}s| calculator {1}: {2}/{3}...".format(time.time()-t0, name, i, NUM))
        # set all multipliers
        gas.set_multiplier(1.0) 
        for j,factor in enumerate(props['factor']):
            gas.set_multiplier(factor, j)

        props['fs'] = get_fs(gas, props)
        # save data
        if props['fs'] != 0:
            data_dict['props'] = deepcopy(props)
            acquireLock(maxt, size, rank)
            json_writer(file_name, data_dict)

if __name__=="__main__":
    # = = = = = = = = = =
    # calculating
    # general conditions
    phi_arr=np.linspace(0.5,1.5,11)
    P_arr=[1,5,20]
    conditions = get_conditions(phi_arr=phi_arr,P_arr=P_arr,T_arr=[300])
    
    # props['width'] = 0.05
    # # for m,name,mech in mechs:
    # m,name,mech = mechs[3]

    # gas = load_mech(mech)
    # uncetainty_factors = load_uncertainty(mech, UF=UF)
    # file_name = comp_dir+name+"_fs_sample.json"
    # # if not checkexists(file_name, delete=False):
    # print("Calculating data of %s on c%d"%(name,rank))

    # props_arr = get_props_arr(props, conditions)
    
    # # random samplings
    # props_arr = []
    # for P in P_arr:
    #     for s in range(100):
    #         phi = np.random.rand()+0.5
    #         props['P'] = P
    #         props['T'] = 300
    #         props['phi'] = phi
    #         props['factor'] = list(get_factors(uncetainty_factors))
    #         props_arr.append(deepcopy(props))

    # l,i = [int(li) for li in np.linspace(0, len(props_arr), num=size+1)],rank
    # calculator('c%d'%i, gas, props_arr[l[i]:l[i+1]], file_name)

    # # else:
    # #     new_props_arr = check_comp(file_name,conditions,delete=False)
    # #     print("Data of %s exists. Need %d calculations."%(name, len(new_props_arr)))
    # #     if len(new_props_arr):
    # #         l,i = [int(li) for li in np.linspace(0, len(new_props_arr), num=size+1)],rank
    # #         calculator('c%d'%i, gas, new_props_arr[l[i]:l[i+1]], file_name)

    # if rank!=0: exit(0)

    # = = = = = = = = = =
    # showing results
    DATA = {}
    maxA ,minA = -1e10, 1e10
    mech_names = ['DME55', 'DME42', 'DME40', 'DME30']

    fig = plt.figure(figsize=c2i(12,9))
    
    for p,P in enumerate(P_arr):
        print("P: %.1f"%P)
        for m,name,mech in mechs:
            print("mech: %s"%name)
            # loading data
            props_arr = [d['props'] for d in json_reader(comp_dir+name+"_fs.json")]
            
            # load data in array T and sorted by 1000/T
            props_parr = [props for props in props_arr if props['P']==P]
            props_parr = sorted(props_parr, key=lambda props:props['phi'])
            phi_arr = [props['phi'] for props in props_parr]
            data_arr = [props['fs']*100 for props in props_parr]
            
            # plot
            plt.plot(phi_arr, data_arr, color_arr[m]+line_arr[m]+symbol_arr[m], fillstyle='none')
            print(*(map('{:5.2f}'.format, phi_arr)), sep=", ")
            print(*(map('{:5.2f}'.format, data_arr)), sep=", ")
            print()
            maxA,minA = max(max(data_arr)+1,maxA),min(min(data_arr)-1,minA)
    
    # = = = = = = = = = =
    # figure setting
    plt.xlabel("Equivalence Ratio")
    plt.ylabel("Laminar Flame Speed, cm/s")
    plt.legend(mech_names)
    plt.ylim([0, 100])
    fig.subplots_adjust(left=0.07,bottom=0.09,top=0.98,right=0.95,hspace=0.,wspace=0.)
    save_figure(fig, path=figs_dir+'compare_fs.png')
    plt.show()
