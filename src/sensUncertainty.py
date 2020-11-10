from utils import *
from check import *

mech = "mech/DME/DMEzhao.cti"
samplingSize = 5
props['T'], props['P'], props['phi'] = 1200., 1., 1.
props['pdiff'] = pdiff

gas = load_mech(mech)
uncetainty_factors = load_uncertainty(mech, UF=5)

print("Original sensitivity")
props['factor'] = list(np.ones(gas.n_reactions))
props['idt'], tdata = get_ign_sens_bf(gas, props)
tdata = normalize(tdata)
rank = np.argsort(-np.abs(tdata))
pos = accum(tdata[rank], sum_limit=0.99)
rank = rank[:15]
x = np.array(range(0,len(rank)))
xtick = np.array(['R%d'%r for r in rank])

print(rank)
v = tdata[rank]
plt.bar(x, tdata[rank], width=1/(samplingSize+1))
#plt.show()

for s in range(1, samplingSize): # random
    print("Sampling: ", s, end=' ')
    props['factor'] = list(get_factors(uncetainty_factors))
    props['idt'], tdata = get_ign_sens_bf(gas,props)
    tdata = normalize(tdata)
    if v[0]*tdata[0]<0: tdata = -tdata
    plt.bar(x+s/(samplingSize+1), tdata[rank], width=1/(samplingSize+1))
    print(np.dot(v, tdata[rank]), np.sqrt(np.sum(tdata[rank]**2)))

plt.ylim([-1,1])
plt.xticks(x+.5, xtick, rotation=45)
plt.savefig("data/DME_sens/sensUncertainty.png", 
    bbox_inches='tight', transparent=True, dpi=300)
plt.show()