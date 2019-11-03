from utils import *
sys.path.append("dep/pymars")
import soln2cti_factor

m,name,mech = mechs[0]

gas = ct.Solution(mech)
# UF = load_uncertainty(mech.split('.')[0]+'.txt')
UF = load_uncertainty(mech, UF=5)

factors =  get_factors(UF)
gas = set_factors(gas, factors)

print(factors)
gas.name = "H2fake"
soln2cti_factor.write(gas)
