import sys
sys.path.append("src")
from utils import *

checkpath(figs_dir)

x = [14,24,38] # DME2000
# x = [33,34,73] # DME

with open('tools/curv_sc2000.txt') as f:
    data = []
    for line in f.readlines():
        lin = line.split(' ')
        lin[-1] = float(lin[-1].split('\n')[0])*0.01

        li = [float(l) for l in lin if l!='']
        data.append(li)
    data = np.transpose(data)

    fig = plt.figure(figsize=(6,4.5))

    ax1 = fig.add_subplot(111)
    ax1.plot(data[0],data[1],'b')
    ax1.set_ylabel(r'$N_{S}$')
    ax1.set_xlabel('threshold '+r'$\varepsilon$')
    ax1.set_title("Mechanism reduction using DRG")

    ax2 = ax1.twinx()  # this is the important function
    ax2.plot(data[0],data[2],'r')
    ax2.plot(data[0][x],data[2][x],'ro')
    ax2.set_yscale(r'log')
    # ax2.set_xlim([0, np.e])
    ax2.set_ylabel(r'$err_{max}$')
    plt.savefig(figs_dir+'DRG.png', format='png', bbox_inches='tight',
                    transparent=True, dpi=600)
    plt.show()
