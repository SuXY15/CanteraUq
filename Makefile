# # Name Rules
# Mechanism Reduction: mr
# Dimension Reduction: dr
# Sensitivity Analysis: sa
# Monte Carlo Sampling: mc

mech = "DME2000"
server = "False"

# configuration
config:
	python3 src/config.py ${mech}

# DRG mechanism reduction
DRG:
	python3 src/DRG.py

# curv drg curv
curv:
	python3 tools/curv.py ${mech}

# IDT mechanism comparation
compare:
ifdef N
	mpiexec -N ${N} python3 src/compare.py ${mech}
else
	python3 src/compare.py ${mech}
endif

# mechanism pdiff calibration
calibration:
ifdef N
	mpiexec -N ${N} python3 src/calibration.py ${mech}
else
	python3 src/calibration.py ${mech}
endif

# mechanism sensitivity analysis
sensAnalyze:
ifdef N
	mpiexec -N ${N} python3 src/sensAnalyze.py ${mech}
else
	python3 src/sensAnalyze.py ${mech}
endif

# mechanism check sensitivity infinty
checkInf:
	python3 src/checkInf.py ${mech}

# active subspace
subspace_sampling:
ifdef N
	mpiexec -N ${N} python3 src/activeSubspace.py ${mech} sampling
else
	python3 src/activeSubspace.py ${mech} sampling
endif
subspace_show:
	python3 src/activeSubspace.py ${mech} show
subspace_generate:
ifdef N
	mpiexec -N ${N} python3 src/activeSubspace.py ${mech} generate
else
	python3 src/activeSubspace.py ${mech} generate
endif

# response surfaces
response_anntrain:
	python3 src/respSurface.py ${mech} anntrain
response_rawtrain:
	python3 src/respSurface.py ${mech} rawtrain
response_train:
	python3 src/respSurface.py ${mech} train
response_predict:
	python3 src/respSurface.py ${mech} predict
response_distribution:
	python3 src/respSurface.py ${mech} distribution

# propagation
type = 3
propagation:
	python3 src/propagation.py ${mech} ${type}

# check
check:
	python3 src/check.py ${mech}

# tools: clean
clean:
	rm data/tmp/*
