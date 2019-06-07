#
# For comparation of detailed mechanism and skeletal mechanism
#
from utils import *
from check import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def Tcurv_mech(props, ax1setting, ax2setting):
    fig = plt.figure("%.0f"%props['T'], figsize=c2i(12,9))
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
        ax1.plot(data['t'], data['T'], color_arr[m]+'-')
        if ax2setting is not None:
            ax2.plot(data['t'], data['T'], color_arr[m]+'-')
    for m,name,mech in mechs:
        gas = load_mech(mech)
        data = getTcurv(gas, props)
        pos,idt = diffMax(data['t'], data['T'])
        ax1.plot(idt, data['T'][pos], color_arr[m]+'o', fillstyle='none')
        legend.append(name)
    ax1.set_xlim(ax1t)
    ax1.set_ylim(ax1T)
    if ax2setting is not None:
        ax2.set_xlim(ax2t); # ax2.set_xticks(ax2t)
        ax2.set_ylim(ax2T); # ax2.set_yticks(ax2T)
        ax1.plot(ax2t,[ax2T[0],ax2T[0]],'k-',lw=0.5); ax1.plot([ax2t[0],ax2t[0]],ax2T,'k-',lw=0.5)
        ax1.plot(ax2t,[ax2T[1],ax2T[1]],'k-',lw=0.5); ax1.plot([ax2t[1],ax2t[1]],ax2T,'k-',lw=0.5)
    ax1.legend(legend, frameon=False, ncol = 2, loc='upper center')
    ax1.set_xlabel(r'Time / s')
    ax1.set_ylabel(r'Temperature / K')
    save_figure(fig, figs_dir+'compare_Tcurv_mech_%.0f.eps'%props['T'])

if __name__=="__main__":
    if mech_arr[0]=="DMEzhao":
        # Low temperature
        props['phi'],props['P'],props['T'] = 1.0, 1.0, 600.0
        ax1t=[0.00,0.15]; ax1T=[500,3500]
        ax2t=[0.09,0.10]; ax2T=[600,1000]
        ax1pos=[0.10, 0.10, 0.80, 0.80]
        ax2pos=[0.25, 0.25, 0.20, 0.35]
        Tcurv_mech(props, [ax1t, ax1T, ax1pos], [ax2t, ax2T, ax2pos])

        # High temperature
        props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200.0
        ax1t=[0.000, 0.002]; ax1T=[1000, 3500]
        ax1pos=[0.10, 0.10, 0.80, 0.80]
        Tcurv_mech(props, [ax1t, ax1T, ax1pos], None)
    if mech_arr[0]=="DME2000":
        # Low temperature
        props['phi'],props['P'],props['T'] = 1.0, 1.0, 600.0
        ax1t=[0.00,0.25]; ax1T=[500,3500]
        ax2t=[0.13,0.16]; ax2T=[600, 900]
        ax1pos=[0.10, 0.10, 0.80, 0.80]
        ax2pos=[0.25, 0.35, 0.35, 0.30]
        Tcurv_mech(props, [ax1t, ax1T, ax1pos], [ax2t, ax2T, ax2pos])

        # High temperature
        props['phi'],props['P'],props['T'] = 1.0, 1.0, 1200.0
        ax1t=[0.000, 0.002]; ax1T=[1000, 3500]
        ax1pos=[0.10, 0.10, 0.80, 0.80]
        Tcurv_mech(props, [ax1t, ax1T, ax1pos], None)

    plt.show()