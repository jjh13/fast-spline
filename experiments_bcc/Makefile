###########################################################
# Cartesian
############################################################
# 16 points in the support, supported on [-2, 2]^3
bcc_tp_lin_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_tp_lin 16 2 > bcc_tp_lin_cpu.experiments.mk

# 16 points in the support, supported on [-1.5, 1.5]^3
bcc_tp_quad_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_tp_quad 54 4 > bcc_tp_quad_cpu.experiments.mk

# 16 points in the support, supported on [-2, 2]^3
bcc_rdodecl_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_rdodecl 4 2 > bcc_rdodecl_cpu.experiments.mk

# 27 points in the support, supported on [-1.5, 1.5]^3
bcc_rdodecq_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_rdodecq 32 4 > bcc_rdodecq_cpu.experiments.mk

# 32 points in the support, supported on [-2, 2]^3
bcc_toctq_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_toctq 30 2 > bcc_toctq_cpu.experiments.mk

# 53 points in the support, supported on ~[-2, 2]^3
bcc_v2_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_v2 8 2 > bcc_v2_cpu.experiments.mk

bcc_v3_cpu.experiments.mk: make_experiment.py
	python make_experiment.py bcc_v3 27 3 > bcc_v3_cpu.experiments.mk

all: bcc_tp_lin_cpu.experiments.mk bcc_tp_quad_cpu.experiments.mk bcc_rdodecl_cpu.experiments.mk bcc_rdodecq_cpu.experiments.mk \
     bcc_toctq_cpu.experiments.mk bcc_v3_cpu.experiments.mk bcc_v2_cpu.experiments.mk

