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
        props['idt'] = get_ign(gas, props)
        # save data
        if props['idt'] != 0:
            data_dict['props'] = deepcopy(props)
            acquireLock(maxt, size, rank)
            json_writer(file_name, data_dict)

if __name__=="__main__":
    # = = = = = = = = = =
    # calculating
    for m,name,mech in mechs:
        gas = load_mech(mech)
        file_name = comp_dir+name+".json"
        # check if exists or not
        if not checkexists(file_name, delete=False):
            print("Calculating data of %s on c%d"%(name,rank))
            props_arr = get_props_arr(props, conditions)
            l,i = [int(li) for li in np.linspace(0, len(props_arr), num=size+1)],rank
            calculator('c%d'%i, gas, props_arr[l[i]:l[i+1]], file_name)
        else:
            new_props_arr = check_comp(name,mech,delete=False)
            print("Data of %s exists. Need %d calculations."%(name, len(new_props_arr)))
            if len(new_props_arr):
                l,i = [int(li) for li in np.linspace(0, len(new_props_arr), num=size+1)],rank
                calculator('c%d'%i, gas, new_props_arr[l[i]:l[i+1]], file_name)

    if rank!=0: exit(0)

    # = = = = = = = = = =
    # showing results
    IDT_DATA = {}
    maxA ,minA = -1e10, 1e10
    maxB ,minB = -1e10, 1e10
    figA, AX = get_sub_plots(num="IDT comparasion")
    figB, BX = get_sub_plots(num="Relative error")

    for m,name,mech in mechs:
        # loading data
        props_arr = [d['props'] for d in json_reader(comp_dir+name+".json")]
        rela_err = []
        for i,phi in enumerate(phi_arr):
            if m==0: IDT_DATA[phi] = {}
            for p,P in enumerate(P_arr):
                # load data in array T and sorted by 1000/T
                props_parr = [props for props in props_arr if props['P']==P and props['phi']==phi]
                props_parr = sorted(props_parr, key=lambda props:1000./props['T'])
                Temp_arr = np.array([props['T'] for props in props_parr if props['T'] in T_arr])
                idt_arr = np.array([props['idt'] for props in props_parr if props['T'] in T_arr])

                # plot IDT
                AX[i,p].plot(1000./Temp_arr,idt_arr,color_arr[m]+'-'+symbol_arr[m])
                maxA,minA = max(max(idt_arr)*32,maxA),min(min(idt_arr)/2,minA)
                
                # error figure
                if m==0: IDT_DATA[phi][P] = idt_arr
                else:
                    diff = np.abs(idt_arr-IDT_DATA[phi][P])/IDT_DATA[phi][P]
                    BX[i,p].plot(1000./Temp_arr,diff,color_arr[m]+'-'+symbol_arr[m])
                    maxB,minB = max(max(diff)*32,maxB),min(min(diff)/2,minB)
                    rela_err += list(diff)
        if len(rela_err):
            cprint("%s, max: %.3f%% mean: %.3f%%"%(name, np.max(rela_err)*100, np.mean(rela_err)*100))

    print(mech_arr)
    # = = = = = = = = = =
    # figure setting
    set_sub_plots(AX, xlabel=r'$1000/T$ (K$^{-1}$)', ylabel=r'$\tau$ (s)',
                    legend=mech_arr,ylim=[minA,maxA], yscale=r'log')
    set_sub_plots(BX, xlabel=r'$1000/T$ (K$^{-1}$)', ylabel=r'$(\tau_s-\tau_d)/\tau_d$',
                    legend=mech_arr[1:],ylim=[minB,maxB], yscale=r'log')
    save_figure(figA, path=figs_dir+'compare_IDT.eps')
    save_figure(figB, path=figs_dir+'compare_err.eps')
    plt.show()
