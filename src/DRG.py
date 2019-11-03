from utils import *

sys.path.append("dep/pymars")
from argparse import ArgumentParser
from os.path import splitext
from warnings import warn

# local imports
import soln2cti
from drg import run_drg
from sensitivity_analysis import run_sa

# ===========================
# setting
model_file = mech_dir+mech_arr[0]+".cti"
conditions = mech_dir+"conditions.txt"
fuel = props['fuel'].split(':')[0]
target_species = [fuel,'O2']
retained_species = [fuel,'O2','N2','H2O']
error = 5.
method = "DRG"
epsilon_star = 0.01

# ===========================
# preparing
solution_object = ct.Solution(model_file)
final_error = [0]

def test():
    thresh = np.array(range(1,50))*0.01
    # running drg
    reduced_model, [errs, nums, dels] = run_drg(solution_object, conditions, error,
            target_species, retained_species, model_file, final_error, thresh)
    with open(mech_dir+'curv.txt','w') as f:
        for i,thr in enumerate(thresh):
            f.writelines("%.3f %4d %8.3f %s\r\n"%(thr, nums[i], errs[i], str(dels[i])))

def reduce(thresh):
    # running drg
    reduced_model, [errs, nums, dels] = run_drg(solution_object, conditions, error,
            target_species, retained_species, model_file, final_error, thresh)
    ori_s = [s.name for s in solution_object.species()]
    new_s = [s.name for s in reduced_model.species()]
    del_s = [s for s in ori_s if s not in new_s]
    cprint("\nFinal counts: sp:%d, reac:%d"%(len(new_s),len(reduced_model.reactions())),'g')
    cprint("Species deleted:"+" ".join(del_s),'g')

    reduced_model.name = mech_dir+"DMEsk%d"%len(new_s)
    output_file = soln2cti.write(reduced_model)

def curv(idx=[]):
    with open(mech_dir+'curv.txt') as f:
        data = []
        for line in f.readlines():
            lin = line.split(' ')
            lin[-1] = lin[-1].split('\n')[0]
            data.append([float(l) for l in lin if l!='' and l[0]>='0' and l[0]<='9' ])
    data = np.array(data)
    data[:,2] /= 100

    fig = plt.figure(figsize=c2i(12,9))

    ax1 = fig.add_subplot(111)
    ax1.plot(data[:,0], data[:,1], 'k-')
    ax1.set_ylabel(r'$N_{S}$')
    ax1.set_xlabel(r'$\varepsilon$')

    ax2 = ax1.twinx()
    ax2.plot(data[:,0], data[:,2], 'r--')
    ax2.plot(data[idx,0],data[idx,2],'ro',ms=4)
    ax2.set_yscale(r'log')
    ax2.set_ylabel(r'$err_{\max}$')
    if len(idx)==0:
        save_figure(fig, path=figs_dir+'DRG_curv.eps')
    else:
        save_figure(fig, path=figs_dir+'DRG.eps')
    plt.show()

if __name__ == "__main__":
    if sys.argv[2] == "test":
        test()
    if sys.argv[2] == "reduce":
        reduce([float(sys.argv[3])])
    if sys.argv[2] == "curv":
        if len(sys.argv)>3:
            idx = [i-1 for i in json.loads(sys.argv[3])]
        else:
            idx = []
        curv(idx)
