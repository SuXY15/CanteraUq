#
# For comparation of detailed mechanism and skeletal mechanism
#
import sys; sys.path.append("src")
from utils import *
from check import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if __name__ == '__main__':    
    m,name,mech = mechs[0]
    gas = load_mech(mech)

    props['T'] = 650.
    fig = plt.figure("property_Tcurv_T",figsize=c2i(12,9))
    for i,props['phi'] in enumerate(phi_arr):
        data = getTcurv(gas, props)
        plt.plot(data['t'], data['T'], color_arr[i]+'-')
    legend = []
    for i,props['phi'] in enumerate(phi_arr):
        data = getTcurv(gas, props)
        pos,idt = diffMax(data['t'], data['T'])
        plt.plot(idt, data['T'][pos], color_arr[i]+'o', fillstyle='none')
        legend.append(r"$\phi = %.1f $"%props['phi'])
    plt.legend(legend, frameon=False)
    plt.xlim([0.015,0.025])
    plt.ylim([600,3200])
    plt.xlabel(r'Time / s')
    plt.ylabel(r'Temperature / K')
    save_figure(fig, figs_dir+'property_Tcurv_T.png')

    fig = plt.figure("property_IDT",figsize=c2i(12,9))
    for p,props['P'] in enumerate(P_arr):
        idt_arr = []
        for props['T'] in T_arr:
            idt_arr.append(get_ign(gas, props))
        plt.plot(1000./np.array(T_arr), idt_arr, color_arr[p]+'-o')
    plt.legend(["P=%.0fatm"%P for P in P_arr], frameon=False)
    plt.yscale(r'log')
    plt.xlabel(r'$1000/T [K^{-1}]$')
    plt.ylabel(r'$IDT [s]$')
    save_figure(fig, figs_dir+'property_IDT.png')

    plt.show()
    