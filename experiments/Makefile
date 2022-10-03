###########################################################
# FCC
############################################################
# 16 points in the support, supported on [-2, 2]^3
fcc_tp_lin_gpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_tp_lin 32 2 --device GPU > fcc_tp_lin_gpu.experiments.mk

# 16 points in the support, supported on [-1.5, 1.5]^3
fcc_tp_quad_gpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_tp_quad 64 4 --device GPU > fcc_tp_quad_gpu.experiments.mk

# 16 points in the support, supported on [-2, 2]^3
fcc_6dir_gpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_6dir 16 2 --device GPU > fcc_6dir_gpu.experiments.mk

# 53 points in the support, supported on ~[-2, 2]^3
fcc_v2_gpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_v2 8 2 --device GPU > fcc_v2_gpu.experiments.mk

fcc_v3_gpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_v3 27 3 --device GPU > fcc_v3_gpu.experiments.mk

gpu: fcc_tp_lin_gpu.experiments.mk  fcc_tp_quad_gpu.experiments.mk fcc_6dir_gpu.experiments.mk fcc_v2_gpu.experiments.mk \
     fcc_v3_gpu.experiments.mk

# 16 points in the support, supported on [-2, 2]^3
fcc_tp_lin_cpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_tp_lin 32 2 > fcc_tp_lin_cpu.experiments.mk

# 16 points in the support, supported on [-1.5, 1.5]^3
fcc_tp_quad_cpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_tp_quad 64 4 > fcc_tp_quad_cpu.experiments.mk

# 16 points in the support, supported on [-2, 2]^3
fcc_6dir_cpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_6dir 16 2 > fcc_6dir_cpu.experiments.mk

# 53 points in the support, supported on ~[-2, 2]^3
fcc_v2_cpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_v2 8 2 > fcc_v2_cpu.experiments.mk

fcc_v3_cpu.experiments.mk: make_experiment.py
	python make_experiment.py fcc_v3 27 3 > fcc_v3_cpu.experiments.mk

cpu: fcc_tp_lin_cpu.experiments.mk  fcc_tp_quad_cpu.experiments.mk fcc_6dir_cpu.experiments.mk fcc_v2_cpu.experiments.mk \
     fcc_v3_cpu.experiments.mk

all: cpu gpu
