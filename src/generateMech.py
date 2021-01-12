from utils import *
# sys.path.append("dep/pymars")
# import soln2cti_factor

print(mechs[0])
m,name,mech = mechs[0]

gas = ct.Solution(mech)
# UF = load_uncertainty(mech.split('.')[0]+'.txt')
UF = load_uncertainty(mech, UF=5)

factors =  get_factors(UF)
gas = set_factors(gas, factors)

for i,factor in enumerate(factors):
    print("R%-3d %10.5f %10.5f"%(i+1, np.log(factor), factor))

UF = load_uncertainty(mech[:-3]+'txt', UF=0.0)
print(np.log(factors)*np.log(UF)*3.)
gas.name = "H2fake"
soln2cti_factor.write(gas)


