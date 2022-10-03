##############################################################
# Path setup
# You should change the SAGEPATH to the local distribution
# of your sage installation. Also, this Makefile will expand
# a good amount of files. If you're strapped for space on
# your local drive, consider setting the cache to a larger
# volume (I have mine symlinked to a path on a large NAS)
##############################################################
SAGEPATH := /mnt/home/john/SageMath
SAGE := $(SAGEPATH)/sage
# CACHE := ncache
CACHE := /mnt/home/john/research/fast_spline_cache
TTRIPLE := x86_64-unknown-linux-gnu
XTRIPLE := aarch64-unknown-linux-gnu
PTXTRIPLE := nvptx64-nvidia-cuda
MCPU := znver2
XCPU := cortex-a72
LANECOUNT := 4
TA := $(TTRIPLE)_$(MCPU)
TX := $(XTRIPLE)_$(XCPU)
CPUFEATS :=
LLC := llc
LLAS := llvm-as
LLC_FLAGS := --enable-no-infs-fp-math \
             --enable-no-nans-fp-math \
             --enable-no-trapping-fp-math \
             --enable-unsafe-fp-math \
             --mcpu=$(MCPU) \
             -relocation-model=pic

NVCC := nvcc
NVSM := sm_75

OPT := opt
OPT_FLAGS := -adce -die -dse -globalopt -instcombine -O3
CFLAGS := -Winline
NVCC_CFLAGS :=
AS := as
GPP := clang++

shell = $*
shell_name = $(shell echo '$*' | tr '[:lower:]' '[:upper:]')
shell_version = ${${shell_name}_VERSION}

FSCMD := fastspline.sage

SPLINE_DIR := $(CACHE)/spline
LATTICE_DIR := $(CACHE)/lattice
RHO_DIR := $(CACHE)/rho
SS_DIR := $(CACHE)/spline_space
LLVM_DIR := $(CACHE)/llvmir
RESULT_DIR := $(CACHE)/results

# Extensions of the outputs of the fastspline.sage script
SSEXT := object.json.ss
SPEXT := object.json.spline
LEXT := object.json.lattice
RHOEXT := object.json.rho
LIREXT := generated_$(TTRIPLE).ll
LOEXT := generated_$(TTRIPLE).opt.ll
ASEXT := generated_$(TTRIPLE).opt.s
GOEXT := generated_$(TTRIPLE).opt.o
BCEXT := generated_$(TTRIPLE).opt.bc
LIREXT := generated_$(TTRIPLE).ll

PLOEXT := generated_$(PTXTRIPLE).opt.ll
PASEXT := generated_$(PTXTRIPLE).opt.s
PGOEXT := generated_$(PTXTRIPLE).opt.o
PBCEXT := generated_$(PTXTRIPLE).opt.bc
PTXEXT := generated_$(PTXTRIPLE).ptx
PLIREXT := generated_$(PTXTRIPLE).ll

# Define all the splines --- each one of these has their own command, so manually
# check the rules below for the specific invocation of the fastspline.sage command
SPLINES := cctp_lin \
           cctp_quad \
           cctp_cubic \
           tp2_lin \
           tp2_quad \
           tp2_cubic \
           rdodec_linear_bcc \
           rdodec_quartic_bcc \
           toct_quintic_bcc \
           fcc6dir \
           zp_element_cc \
           v_fcc2 \
           v_fcc3 \
           v_fcc4 \
           v_bcc2 \
           v_bcc3 \
           v_bcc4

RHOS := rho_cctp0 \
        rho_cctp1 \
        rho_bcc_ind \
        rho_bcc0 \
        rho_fcc0

# Define all the lattices --- each one of these has their own command, so manually
# check the rules below for the specific invocation of the fastspline.sage command
LATTICES :=  cc cp qc fcc bcc d4 d4s

# Define all the lattice-basis combinations we're investigating in this work
# The convention for these is:
#    spline_space_name:spline_name:lattice:rho
# spline space name must start with the lattice, i.e. cc_rest_of_name
SPLINE_SPACES := cc_tplin:cctp_lin:cc:rho_cctp0 \
                 cc_tpquad:cctp_quad:cc:rho_cctp1 \
                 cc_tpcubic:cctp_cubic:cc:rho_cctp0 \
                 cc_zp_element:zp_element_cc:cc:rho_cctp1 \
                 cc_6dir:fcc6dir:cc:rho_cctp1 \
                 cc_rdodecl:rdodec_linear_bcc:cc:rho_cctp1 \
                 cc_rdodecq:rdodec_quartic_bcc:cc:rho_cctp1 \
                 cc_toctq:toct_quintic_bcc:cc:rho_cctp1 \
                 cc_fcc_v2:v_fcc2:cc:rho_cctp1 \
                 cc_fcc_v3:v_fcc3:cc:rho_cctp1 \
                 cc_bcc_v2:v_bcc2:cc:rho_cctp1 \
                 cc_bcc_v3:v_bcc3:cc:rho_cctp1 \
                 bcc_tp_lin:tp2_lin:bcc:rho_bcc_ind \
                 bcc_tp_quad:tp2_quad:bcc:rho_bcc_ind \
                 bcc_rdodecl:rdodec_linear_bcc:bcc:rho_bcc_ind \
                 bcc_rdodecq:rdodec_quartic_bcc:bcc:rho_bcc_ind \
                 bcc_toctq:toct_quintic_bcc:bcc:rho_bcc_ind \
                 bcc_v2:v_bcc2:bcc:rho_bcc_ind \
                 bcc_v3:v_bcc3:bcc:rho_bcc_ind \
                 fcc_tp_lin:tp2_lin:fcc:rho_fcc1 \
                 fcc_tp_quad:tp2_quad:fcc:rho_fcc1 \
                 fcc_v2:v_fcc2:fcc:rho_fcc1 \
                 fcc_v3:v_fcc3:fcc:rho_fcc1 \
                 fcc_6dir:fcc6dir:fcc:rho_fcc1

# bcc_lin_rdodec:rdodec_linear_bcc:bcc:rho_bcc0 \
#                  fcc_v2:v_fcc2:fcc:rho_fcc0 \
#                  fcc_v3:v_fcc3:fcc:rho_fcc1


#                 cc_tplin:cctp_lin:cc:rho_cctp0 \
#                  cc_tpquad:cctp_quad:cc:rho_cctp1 \
#                  cc_tpcubic:cctp_cubic:cc:rho_cctp0 \
#                  bcc_lin_rdodec:rdodec_linear_bcc:bcc:rho_bcc0
#                  bcc_quart_rdodec:rdodec_quartic_bcc:bcc:rho_bcc0 \
#                  bcc_quint_toct:toct_quintic_bcc:bcc:rho_bcc0 \
#                  fcc_6dir:fcc6dir:fcc:rho_fcc0 \
#                  cc_6dir:fcc6dir:cc:rho_cctp1 \
#                  bcc_v2:v_bcc2:bcc:rho_bcc0 \
#                  fcc_v2:v_fcc2:fcc:rho_fcc0
#                  fcc_v3:v_fcc3:fcc:rho_fcc0 \
#                  bcc_v3:v_bcc3:bcc:rho_bcc0 \
#                  fcc_v4:v_fcc4:fcc:rho_fcc0
#                  bcc_v4:v_bcc4:bcc:rho_bcc0

##############################################################
# Rule start
##############################################################
all: derive_splines derive_lattices derive_ss derive_rho
derive_splines: $(foreach S, $(SPLINES), $(SPLINE_DIR)/$(S)/$(SPEXT))
derive_lattices: $(foreach L, $(LATTICES), $(LATTICE_DIR)/$(L)/$(LEXT))
derive_rho: $(foreach R, $(RHOS), $(RHO_DIR)/$(R)/$(RHOEXT))

clean_splines:
	rm -rf $(SPLINE_DIR)

clean_lattice:
	rm -rf $(SPLINE_DIR)

clean_spline_space:
	rm -rf $(SS_DIR)

clean_llvm:
	rm -rf $(LLVM_DIR)

clean: clean_splines clean_lattice clean_llvm clean_spline_space

##############################################################
# Functionality testing
##############################################################
test: test/test_scalar test/test_vec4 test/baked_test test/baked_test_vec

test/test_scalar: $(SS_DIR)/cc_tplin/$(SSEXT) instrumentation/drivers/abstract_scalar.c
	$(SAGE) $(FSCMD) codegen \
	    -ss $(SS_DIR)/cc_tplin/$(SSEXT) \
	    -o $(LLVM_DIR)/test/test_scalar/generated.ll
	$(OPT) $(LLVM_DIR)/test/test_scalar/generated.ll $(OPT_FLAGS) -S -o $(LLVM_DIR)/test/test_scalar/generated.opt.ll
	$(LLC) $(LLVM_DIR)/test/test_scalar/generated.opt.ll $(LLC_FLAGS) \
	       -o $(LLVM_DIR)/test/test_scalar/generated.opt.s
	$(AS) $(LLVM_DIR)/test/test_scalar/generated.opt.s -o $(LLVM_DIR)/test/test_scalar/generated.opt.o
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_CC -DDIMENSION=3 \
	                    instrumentation/drivers/abstract_scalar.c \
	                    $(LLVM_DIR)/test/test_scalar/generated.opt.o \
	                    -o test/test_scalar

test/test_vec4: $(SS_DIR)/cc_tplin/$(SSEXT) instrumentation/drivers/abstract_vec4.c
	$(SAGE) $(FSCMD) codegen \
	    -ss $(SS_DIR)/cc_tplin/$(SSEXT) \
	    -o $(LLVM_DIR)/test/test_vec4/generated.ll \
	    --cg-lanes 4
	$(OPT) $(LLVM_DIR)/test/test_vec4/generated.ll $(OPT_FLAGS) -S -o $(LLVM_DIR)/test/test_vec4/generated.opt.ll
	$(LLC) $(LLVM_DIR)/test/test_vec4/generated.opt.ll $(LLC_FLAGS) \
	       -o $(LLVM_DIR)/test/test_vec4/generated.opt.s
	$(AS) $(LLVM_DIR)/test/test_vec4/generated.opt.s -o $(LLVM_DIR)/test/test_vec4/generated.opt.o
	$(GPP) -march=znver2 -O3 -DLATTICE_CC -DDIMENSION=3 \
	                    instrumentation/drivers/abstract_vec4.c \
	                    $(LLVM_DIR)/test/test_vec4/generated.opt.o \
	                    -o test/test_vec4

test/baked_test: $(SS_DIR)/cc_tplin/$(SSEXT) instrumentation/drivers/coset_raw_scalar.c codegen.sage
	$(SAGE) $(FSCMD) codegen \
	    -ss $(SS_DIR)/cc_tplin/$(SSEXT) \
	    -o $(LLVM_DIR)/test/baked_test/generated.ll \
	    --cg-lanes 1 \
	    --cg-fetches COSET_RAW
	$(OPT) $(LLVM_DIR)/test/baked_test/generated.ll $(OPT_FLAGS) -S -o $(LLVM_DIR)/test/baked_test/generated.opt.ll
	$(LLC) $(LLVM_DIR)/test/baked_test/generated.opt.ll $(LLC_FLAGS) \
	       -o $(LLVM_DIR)/test/baked_test/generated.opt.s
	$(AS) $(LLVM_DIR)/test/baked_test/generated.opt.s -o $(LLVM_DIR)/test/baked_test/generated.opt.o
	$(GPP) -march=znver2 -O3 -DLATTICE_CC -DDIMENSION=3 \
	                    instrumentation/drivers/coset_raw_scalar.c \
	                    $(LLVM_DIR)/test/baked_test/generated.opt.o \
	                    -o test/baked_test


test/baked_test_vec: $(SS_DIR)/cc_tplin/$(SSEXT) instrumentation/drivers/coset_raw_vec.c codegen.sage
	$(SAGE) $(FSCMD) codegen \
	    -ss $(SS_DIR)/cc_tplin/$(SSEXT) \
	    -o $(LLVM_DIR)/test/baked_test_vec/generated.ll \
	    --cg-lanes 4 \
	    --cg-pipeline-depth 1 \
	    --cg-group-size 8 \
	    --cg-fetches COSET_RAW
	$(OPT) $(LLVM_DIR)/test/baked_test_vec/generated.ll $(OPT_FLAGS) -S -o $(LLVM_DIR)/test/baked_test_vec/generated.opt.ll
	$(LLC) $(LLVM_DIR)/test/baked_test_vec/generated.opt.ll $(LLC_FLAGS) \
	       -o $(LLVM_DIR)/test/baked_test_vec/generated.opt.s
	$(AS) $(LLVM_DIR)/test/baked_test_vec/generated.opt.s -o $(LLVM_DIR)/test/baked_test_vec/generated.opt.o
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_CC -DDIMENSION=3 \
	                    instrumentation/drivers/coset_raw_vec.c \
	                    $(LLVM_DIR)/test/baked_test_vec/generated.opt.o \
	                    -o test/baked_test_vec

##############################################################
# Derive all the \rho s that we use in our experiments
##############################################################
$(RHO_DIR)/rho_cctp0/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt indicator --rho-shift "[1/2, 1/2, 1/2]" -o $@

$(RHO_DIR)/rho_cctp1/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt indicator --rho-shift "[0, 0, 0]" -o $@

$(RHO_DIR)/rho_bcc_ind/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt indicator --rho-shift "[0, 0, 0]" -o $@

$(RHO_DIR)/rho_bcc0/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt pp -rpp "[-1,1,1;1,-1,1;1,1,-1]" --rho-shift "[0, 0, 0]" -o $@

$(RHO_DIR)/rho_fcc0/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt indicator --rho-shift "[0,0,1/2]" -o $@

$(RHO_DIR)/rho_fcc1/object.json.rho:
	$(SAGE) $(FSCMD) rho -rt indicator  --rho-shift "[1/2,1/2,1/2]" -o $@

##############################################################
# Derive all the splines that we use in our experiments
##############################################################
$(SPLINE_DIR)/cctp_lin/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[1,0,0,1,0,0;0,1,0,0,1,0;0,0,1,0,0,1]" -o $@

$(SPLINE_DIR)/tp2_lin/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[2,0,0,2,0,0;0,2,0,0,2,0;0,0,2,0,0,2]" -o $@

$(SPLINE_DIR)/cctp_quad/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[1,0,0,1,0,0,1,0,0;0,1,0,0,1,0,0,1,0;0,0,1,0,0,1,0,0,1]" -o $@

$(SPLINE_DIR)/tp2_quad/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[2,0,0,2,0,0,2,0,0;0,2,0,0,2,0,0,2,0;0,0,2,0,0,2,0,0,2]" -o $@

$(SPLINE_DIR)/cctp_cubic/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[1,0,0,1,0,0,1,0,0,1,0,0;0,1,0,0,1,0,0,1,0,0,1,0;0,0,1,0,0,1,0,0,1,0,0,1]" -o $@

$(SPLINE_DIR)/tp2_cubic/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[2,0,0,2,0,0,2,0,0,2,0,0;0,2,0,0,2,0,0,2,0,0,2,0;0,0,2,0,0,2,0,0,2,0,0,2]" -o $@

$(SPLINE_DIR)/rdodec_linear_bcc/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[-1,1,1,-1;1,-1,1,-1;1,1,-1,-1]" -o $@

$(SPLINE_DIR)/rdodec_quartic_bcc/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[-1,1,1,-1, -1,1,1,-1;1,-1,1,-1,1,-1,1,-1;1,1,-1,-1,1,1,-1,-1]" -o $@

$(SPLINE_DIR)/toct_quintic_bcc/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[-1,1,1,-1, 2,0,0;1,-1,1,-1,0,2,0;1,1,-1,-1,0,0,2]" -o $@

$(SPLINE_DIR)/fcc6dir/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[1,-1,1,1,0,0;1,1,0,0,1,-1;0,0,1,-1,1,1]" -o $@

$(SPLINE_DIR)/zp_element_cc/object.json.spline:
	$(SAGE) $(FSCMD) derive -sm "[-1,1,1,-1, 1,0,0;1,-1,1,-1,0,1,0;1,1,-1,-1,0,0,1]" -o $@

$(SPLINE_DIR)/v_fcc2/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s fcc2 -o $@

$(SPLINE_DIR)/v_fcc3/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s fcc3 -o $@

$(SPLINE_DIR)/v_fcc4/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s fcc4 -o $@

$(SPLINE_DIR)/v_bcc2/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s bcc2 -o $@

$(SPLINE_DIR)/v_bcc3/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s bcc3 -o $@

$(SPLINE_DIR)/v_bcc4/object.json.spline:
	$(SAGE) $(FSCMD) derive_voronoi -s bcc4 -o $@

##############################################################
# Generate lattice structures for our experiments
##############################################################
$(LATTICE_DIR)/bcc/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[-1,1,1;1,-1,1;1,1,-1]" -o $@

$(LATTICE_DIR)/cc/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[1,0,0;0,1,0;0,0,1]" -o $@

$(LATTICE_DIR)/cp/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[1,0;0,1]" -o $@

$(LATTICE_DIR)/qc/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[1,1;1,-1]" -o $@

$(LATTICE_DIR)/fcc/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[0,1,1;1,0,1;1,1,0]" -o $@

$(LATTICE_DIR)/d4/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[-1,0,0,0;1,0,-1,-1;0,-1,0,1;0,-1,1,0]" -o $@

$(LATTICE_DIR)/d4s/object.json.lattice:
	$(SAGE) $(FSCMD) define_lattice -lm "[2,0,0,1; 0,2,0,1; 0,0,2,1; 0,0,0,1]" -o $@

##############################################################
# Setup all the Spline Space rules
##############################################################
define SPLINESPACE_RULE
$(SS_DIR)/$(1)/$(SSEXT): $(SPLINE_DIR)/$(2)/$(SPEXT) $(LATTICE_DIR)/$(3)/$(LEXT) $(RHO_DIR)/$(4)/$(RHOEXT)
	$(SAGE) $(FSCMD) compute_spline_space -s $(SPLINE_DIR)/$(2)/$(SPEXT) \
	                 -l $(LATTICE_DIR)/$(3)/$(LEXT) \
	                 -o $(SS_DIR)/$(1)/$(SSEXT) \
	                 -r $(RHO_DIR)/$(4)/$(RHOEXT) \
	                 --ss-disable-reflective False \
	                 --region-cache-path $(SS_DIR)/$(1)/
endef

$(foreach SS,$(SPLINE_SPACES), \
    $(eval SSP = $(word 1,$(subst :, ,$(SS)))) \
    $(eval SP = $(word 2,$(subst :, ,$(SS)))) \
    $(eval L = $(word 3,$(subst :, ,$(SS)))) \
    $(eval R = $(word 4,$(subst :, ,$(SS)))) \
    $(eval $(call SPLINESPACE_RULE,$(SSP),$(SP),$(L),$(R))) \
)

derive_ss: $(foreach SS,$(SPLINE_SPACES), $(eval SSP = $(word 1,$(subst :, ,$(SS))))  $(SS_DIR)/$(SSP)/$(SSEXT))
##############################################################
# Setup all the experiment rules
##############################################################

GPU_EXP :=
GPU_EXP_RES :=

GPU_REN_EXP :=
GPU_REN_RES :=

CPU_EXP := 
CPU_EXP_RES :=
XCPU_EXP :=
XCPU_EXP_RES :=

CPU_VAL :=
CPU_VAL_RES :=
XCPU_VAL :=
XCPU_VAL_RES :=


define SS_GPU_EXPERIMENT

$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLIREXT): $(SS_DIR)/$(1)/$(SSEXT)
	$(SAGE) $(FSCMD) codegen \
        -ss $(SS_DIR)/$(1)/$(SSEXT) \
        -o $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/ \
        --cg-lanes 1 \
        --cg-group-size $(2) \
        --cg-pipeline-depth $(3) \
        --cg-fetches TEX \
        --cg-target nvptx64-nvidia-cuda \
        --cg-func-name reconstruct_$(word 1,$(subst _, ,$(1))) \
        --region-cache-path $(SS_DIR)/$(1)/


$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLOEXT): $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLOEXT)

$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PTXEXT): $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLOEXT)
	$(LLC) -mcpu=$(NVSM) $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PLOEXT) -o $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PTXEXT)

$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o: $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PTXEXT) instrumentation/benchmark/kernel_dummy_fcc.cu
	python asptx.py  \
		$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/$(PTXEXT) \
		instrumentation/benchmark/kernel_dummy_$(word 1,$(subst _, ,$(1))).cu \
		$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o \
		-arch=$(NVSM)

$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda: $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o instrumentation/benchmark/cuda_benchmark.cu
	$(NVCC) instrumentation/benchmark/cuda_benchmark.cu \
	    -x cu \
	    -dc \
		-O3 $(NVCC_CFLAGS) \
		-DLATTICE_$(word 1,$(subst _, ,$(1))) \
		-DDIMENSION=3 \
		-o $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.o \
		-arch=$(NVSM)

	$(NVCC) $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.o \
	    -arch=$(NVSM) \
	    $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o \
	    -o $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda

$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.csv: $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda
	$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda > $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.csv

$(RESULT_DIR)/$1/e_gpu_$(2)_$(3)/benchmark_cuda.csv: $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.csv
	mkdir -p $(RESULT_DIR)/$1/e_gpu_$(2)_$(3)/
	cp $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda.csv $(RESULT_DIR)/$1/e_gpu_$(2)_$(3)/benchmark_cuda.csv

GPU_EXP += $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/benchmark_cuda
GPU_EXP_RES += $(RESULT_DIR)/$1/e_gpu_$(2)_$(3)/benchmark_cuda.csv

endef



define SS_GPU_LIN_EXPERIMENT

$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLIREXT): $(SS_DIR)/$(1)/$(SSEXT)
	$(SAGE) $(FSCMD) codegen \
        -ss $(SS_DIR)/$(1)/$(SSEXT) \
        -o $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/ \
        --cg-lanes 1 \
        --cg-pipeline-depth $(2) \
        --cg-fetches TEX_LIN \
        --cg-target nvptx64-nvidia-cuda \
        --cg-func-name reconstruct_$(word 1,$(subst _, ,$(1))) \
        --region-cache-path $(SS_DIR)/$(1)/


$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLOEXT): $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLOEXT)

$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PTXEXT): $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLOEXT)
	$(LLC) -mcpu=$(NVSM) $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PLOEXT) -o $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PTXEXT)

$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o: $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PTXEXT) instrumentation/benchmark/kernel_dummy_fcc.cu
	python asptx.py  \
		$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/$(PTXEXT) \
		instrumentation/benchmark/kernel_dummy_$(word 1,$(subst _, ,$(1))).cu \
		$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o \
		-arch=$(NVSM)

$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda: $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o instrumentation/benchmark/cuda_benchmark.cu
	$(NVCC) instrumentation/benchmark/cuda_benchmark.cu \
	    -x cu \
	    -dc \
		-O3 $(NVCC_CFLAGS) \
		-DLATTICE_$(word 1,$(subst _, ,$(1))) \
		-DDIMENSION=3 \
		-DTEST_LINEAR_FETCH \
		-o $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.o \
		-arch=$(NVSM)

	$(NVCC) $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.o \
	    -arch=$(NVSM) \
	    -DTEST_LINEAR_FETCH \
	    $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o \
	    -o $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda

$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.csv: $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda
	$(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda > $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.csv

$(RESULT_DIR)/$1/e_gpu_lin_$(2)/benchmark_lin_cuda.csv: $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.csv
	mkdir -p $(RESULT_DIR)/$1/e_gpu_lin_$(2)/
	cp $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda.csv $(RESULT_DIR)/$1/e_gpu_lin_$(2)/benchmark_lin_cuda.csv

GPU_LIN_EXP += $(LLVM_DIR)/$(1)/e_gpu_lin_$(2)/benchmark_lin_cuda
GPU_LIN_EXP_RES += $(RESULT_DIR)/$1/e_gpu_lin_$(2)/benchmark_lin_cuda.csv

endef

define SS_GPU_RENDER

$(LLVM_DIR)/$(1)/vol_render_$(2)_$(3): $(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o instrumentation/application/volren/volumeRender.cpp instrumentation/application/volren/volumeRender_kernel.cu
	$(NVCC) instrumentation/application/volren/volumeRender.cpp \
	    -x cu \
	    -dc \
		-O3 $(NVCC_CFLAGS) \
		-DLATTICE_$(word 1,$(subst _, ,$(1))) \
		-DDIMENSION=3 \
		-I./instrumentation/include \
		-o $(LLVM_DIR)/$(1)/vr_$(2)_$(3).o \
		-arch=$(NVSM)

	$(NVCC) instrumentation/application/volren/volumeRender_kernel.cu \
	    -x cu \
	    -dc \
		-O3 $(NVCC_CFLAGS) \
		-I./instrumentation/include \
		-DLATTICE_$(word 1,$(subst _, ,$(1))) \
		-DDIMENSION=3 \
		-o $(LLVM_DIR)/$(1)/vr_k_$(2)_$(3).o \
		-arch=$(NVSM)

	$(NVCC) -O3 $(NVCC_CFLAGS) \
		-DLATTICE_$(word 1,$(subst _, ,$(1))) \
		-DDIMENSION=3 \
		-arch=$(NVSM) \
		-o $(LLVM_DIR)/$(1)/vol_render_$(2)_$(3) \
		-lGL -lglut \
		$(LLVM_DIR)/$(1)/vr_k_$(2)_$(3).o \
		$(LLVM_DIR)/$(1)/vr_$(2)_$(3).o \
		$(LLVM_DIR)/$(1)/e_gpu_$(2)_$(3)/kernel_$(PTXTRIPLE)_for_$(TTRIPLE).o

$(RESULT_DIR)/$1/vol_render_$(2)_$(3).csv: $(LLVM_DIR)/$(1)/vol_render_$(2)_$(3)
	mkdir -p $(RESULT_DIR)/$1/
	$(LLVM_DIR)/$(1)/vol_render_$(2)_$(3) 2> $(RESULT_DIR)/$1/vol_render_$(2)_$(3).csv

GPU_REN_EXP += $(LLVM_DIR)/$(1)/vol_render_$(2)_$(3)
GPU_REN_RES += $(RESULT_DIR)/$1/vol_render_$(2)_$(3).csv
endef


define SS_CPU_EXPERIMENT
# codegen.sage

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LIREXT): $(SS_DIR)/$(1)/$(SSEXT)
	$(SAGE) $(FSCMD) codegen \
	-ss $(SS_DIR)/$(1)/$(SSEXT) \
	-o $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/ \
	--cg-lanes 1 \
	--cg-group-size $(2) \
	--cg-pipeline-depth $(3) \
	--cg-fetches COSET_RAW \
	--region-cache-path $(SS_DIR)/$(1)/

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LOEXT): $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LOEXT)

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(BCEXT): $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LOEXT)
	$(LLAS) $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LIREXT) \
		-o $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(BCEXT)

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(ASEXT): $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LOEXT)
	$(LLC) $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(LOEXT) $(LLC_FLAGS) \
		-o $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(ASEXT)

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(GOEXT): $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(ASEXT)
	$(AS) $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(ASEXT) \
		-o $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(GOEXT)

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/benchmark_$(TA): $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(GOEXT) $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(BCEXT) instrumentation/benchmark/cpu_scalar_benchmark.cpp
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_scalar_benchmark.cpp \
	                    $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(GOEXT) \
			    -o $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA) ; \
	mkdir -p $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/

$(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/benchmark_$(TX): instrumentation/benchmark/cpu_scalar_benchmark.cpp
	$(GPP) -mcpu=$(XCPU) --target=$(XTRIPLE) -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_scalar_benchmark.cpp \
	                    $(LLVM_DIR)/$(1)/e_scalar_$(2)_$(3)/$(BCEXT) \
			    -o $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX) ; \
	mkdir -p $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/

$(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/: $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA)
	mkdir -p $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/

$(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA).csv: $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/ $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA)
	mkdir -p $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/
	$(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA) > $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA).csv

$(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX).csv: $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/ $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX)
	mkdir -p $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/
	$(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX) > $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX).csv

CPU_EXP += $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA)
CPU_EXP_RES += $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TA).csv
XCPU_EXP += $(LLVM_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX)
XCPU_EXP_RES += $(RESULT_DIR)/$1/e_scalar_$(2)_$(3)/benchmark_$(TX).csv
endef

define SS_CPU_VALIDATION
$(LLVM_DIR)/$(1)/e_val/$(LIREXT): $(SS_DIR)/$(1)/$(SSEXT)
	$(SAGE) $(FSCMD) codegen \
	-ss $(SS_DIR)/$(1)/$(SSEXT) \
	-o $(LLVM_DIR)/$(1)/e_val/ \
	--cg-lanes 1 \
	--cg-group-size 1 \
	--cg-pipeline-depth 1 \
	--cg-fetches COSET_RAW \
	--region-cache-path $(SS_DIR)/$(1)/

$(LLVM_DIR)/$(1)/e_val/$(LOEXT): $(LLVM_DIR)/$(1)/e_val/$(LIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_val/$(LIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_val/$(LOEXT)

$(LLVM_DIR)/$(1)/e_val/$(BCEXT): $(LLVM_DIR)/$(1)/e_val/$(LOEXT)
	$(LLAS) $(LLVM_DIR)/$(1)/e_val/$(LIREXT) \
		-o $(LLVM_DIR)/$(1)/e_val/$(BCEXT)

$(LLVM_DIR)/$(1)/e_val/$(ASEXT): $(LLVM_DIR)/$(1)/e_val/$(LOEXT)
	$(LLC) $(LLVM_DIR)/$(1)/e_val/$(LOEXT) $(LLC_FLAGS) \
		-o $(LLVM_DIR)/$(1)/e_val/$(ASEXT)

$(LLVM_DIR)/$(1)/e_val/$(GOEXT): $(LLVM_DIR)/$(1)/e_val/$(ASEXT)
	$(AS) $(LLVM_DIR)/$(1)/e_val/$(ASEXT) \
		-o $(LLVM_DIR)/$(1)/e_val/$(GOEXT)

$(LLVM_DIR)/$(1)/e_val/testbench_$(TA): $(LLVM_DIR)/$(1)/e_val/$(GOEXT) $(LLVM_DIR)/$(1)/e_val/$(BCEXT) instrumentation/benchmark/cpu_scalar_testbench.cpp
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_scalar_testbench.cpp \
	                    $(LLVM_DIR)/$(1)/e_val/$(GOEXT) \
			    -o $(LLVM_DIR)/$1/e_val/benchmark_$(TA) ; \
	mkdir -p $(RESULT_DIR)/$1/e_val/

$(LLVM_DIR)/$(1)/e_val/testbench_$(TX): instrumentation/benchmark/cpu_scalar_testbench.cpp
	$(GPP) -mcpu=$(XCPU) --target=$(XTRIPLE) -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_scalar_testbench.cpp \
	                    $(LLVM_DIR)/$(1)/e_val/$(BCEXT) \
			    -o $(LLVM_DIR)/$1/e_val/benchmark_$(TX) ; \
	mkdir -p $(RESULT_DIR)/$1/e_val/

$(RESULT_DIR)/$1/e_val/: $(LLVM_DIR)/$1/e_val/testbench_$(TA)
	mkdir -p $(RESULT_DIR)/$1/e_val/


$(RESULT_DIR)/$1/e_val/test_$(TA).gif: $(RESULT_DIR)/$1/e_val/ $(LLVM_DIR)/$1/e_val/testbench_$(TA)
	$(LLVM_DIR)/$1/e_val/benchmark_$(TA) $(2) > $(RESULT_DIR)/$1/e_val/test_data.eval_py
	python render_test.py $(RESULT_DIR)/$1/e_val/test_data.eval_py $(RESULT_DIR)/$1/e_val/test_$(TA).gif

$(RESULT_DIR)/$1/e_val/test_$(TX).gif: $(RESULT_DIR)/$1/e_val/ $(LLVM_DIR)/$1/e_val/testbench_$(TX)
	$(LLVM_DIR)/$1/e_val/benchmark_$(TX) $(2) > $(RESULT_DIR)/$1/e_val/test_data.eval_py
	python render_test.py $(RESULT_DIR)/$1/e_val/test_data.eval_py $(RESULT_DIR)/$1/e_val/test_$(TX).gif

CPU_VAL += $(LLVM_DIR)/$1/e_val/testbench_$(TA)
CPU_VAL_RES += $(RESULT_DIR)/$1/e_val/test_$(TA).gif
XCPU_VAL += $(LLVM_DIR)/$1/e_val/testbench_$(TX)
XCPU_VAL_RES += $(RESULT_DIR)/$1/e_val/test_$(TX).gif
endef

define SS_CPU_VEC_VALIDATION
$(LLVM_DIR)/$(1)/e_val_vec/$(LIREXT): $(SS_DIR)/$(1)/$(SSEXT) codegen.sage
	$(SAGE) $(FSCMD) codegen \
	-ss $(SS_DIR)/$(1)/$(SSEXT) \
	-o $(LLVM_DIR)/$(1)/e_val_vec/ \
	--cg-lanes 4 \
	--cg-group-size 4 \
	--cg-pipeline-depth 1 \
	--cg-fetches COSET_RAW \
	--region-cache-path $(SS_DIR)/$(1)/

$(LLVM_DIR)/$(1)/e_val_vec/$(LOEXT): $(LLVM_DIR)/$(1)/e_val_vec/$(LIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_val_vec/$(LIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_val_vec/$(LOEXT)

$(LLVM_DIR)/$(1)/e_val_vec/$(BCEXT): $(LLVM_DIR)/$(1)/e_val_vec/$(LOEXT)
	$(LLAS) $(LLVM_DIR)/$(1)/e_val_vec/$(LIREXT) \
		-o $(LLVM_DIR)/$(1)/e_val_vec/$(BCEXT)

$(LLVM_DIR)/$(1)/e_val_vec/$(ASEXT): $(LLVM_DIR)/$(1)/e_val_vec/$(LOEXT)
	$(LLC) $(LLVM_DIR)/$(1)/e_val_vec/$(LOEXT) $(LLC_FLAGS) \
		-o $(LLVM_DIR)/$(1)/e_val_vec/$(ASEXT)

$(LLVM_DIR)/$(1)/e_val_vec/$(GOEXT): $(LLVM_DIR)/$(1)/e_val_vec/$(ASEXT)
	$(AS) $(LLVM_DIR)/$(1)/e_val_vec/$(ASEXT) \
		-o $(LLVM_DIR)/$(1)/e_val_vec/$(GOEXT)

$(LLVM_DIR)/$(1)/e_val_vec/testbench_$(TA): $(LLVM_DIR)/$(1)/e_val_vec/$(GOEXT) $(LLVM_DIR)/$(1)/e_val_vec/$(BCEXT) instrumentation/benchmark/cpu_vector_testbench.cpp
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_vector_testbench.cpp \
	                    $(LLVM_DIR)/$(1)/e_val_vec/$(BCEXT) \
			    -o $(LLVM_DIR)/$1/e_val_vec/benchmark_$(TA) ; \
	mkdir -p $(RESULT_DIR)/$1/e_val_vec/

$(LLVM_DIR)/$(1)/e_val_vec/testbench_$(TX): instrumentation/benchmark/cpu_vector_testbench.cpp
	$(GPP) -mcpu=$(XCPU) --target=$(XTRIPLE) -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
	                    -DDIMENSION=3 \
	                    instrumentation/benchmark/cpu_vector_testbench.cpp \
	                    $(LLVM_DIR)/$(1)/e_val_vec/$(BCEXT) \
			    -o $(LLVM_DIR)/$1/e_val_vec/benchmark_$(TX) ; \
	mkdir -p $(RESULT_DIR)/$1/e_val_vec/

$(RESULT_DIR)/$1/e_val_vec/: $(LLVM_DIR)/$1/e_val_vec/testbench_$(TA)
	mkdir -p $(RESULT_DIR)/$1/e_val_vec/

$(RESULT_DIR)/$1/e_val_vec/test_$(TA).gif: $(RESULT_DIR)/$1/e_val_vec/ $(LLVM_DIR)/$1/e_val_vec/testbench_$(TA)
	$(LLVM_DIR)/$1/e_val_vec/benchmark_$(TA) $(2) > $(RESULT_DIR)/$1/e_val_vec/test_data.eval_py
	python render_test.py $(RESULT_DIR)/$1/e_val_vec/test_data.eval_py $(RESULT_DIR)/$1/e_val_vec/test_$(TA).gif

$(RESULT_DIR)/$1/e_val_vec/test_$(TX).gif: $(RESULT_DIR)/$1/e_val_vec/ $(LLVM_DIR)/$1/e_val_vec/testbench_$(TX)
	$(LLVM_DIR)/$1/e_val_vec/benchmark_$(TX) $(2) > $(RESULT_DIR)/$1/e_val_vec/test_data.eval_py
	python render_test.py $(RESULT_DIR)/$1/e_val_vec/test_data.eval_py $(RESULT_DIR)/$1/e_val_vec/test_$(TX).gif

CPU_VAL += $(LLVM_DIR)/$1/e_val_vec/testbench_$(TA)
CPU_VAL_RES += $(RESULT_DIR)/$1/e_val_vec/test_$(TA).gif
XCPU_VAL += $(LLVM_DIR)/$1/e_val_vec/testbench_$(TX)
XCPU_VAL_RES += $(RESULT_DIR)/$1/e_val_vec/test_$(TX).gif
endef

define SS_CPU_VEC_EXPERIMENT
#codegen.sage

$(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LIREXT): $(SS_DIR)/$(1)/$(SSEXT)
	$(SAGE) $(FSCMD) codegen \
	-ss $(SS_DIR)/$(1)/$(SSEXT) \
	-o $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/ \
	--cg-lanes $(LANECOUNT) \
	--cg-group-size $(2) \
	--cg-pipeline-depth $(3) \
	--cg-fetches COSET_RAW \
	--region-cache-path $(SS_DIR)/$(1)/

$(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT): $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LIREXT)
	$(OPT) $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LIREXT) $(OPT_FLAGS) \
		-S -o $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT)

$(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(ASEXT): $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT)
	$(LLC) $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT) $(LLC_FLAGS) \
		-o $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(ASEXT)

$(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(GOEXT): $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(ASEXT)
	$(AS) $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(ASEXT) \
		-o $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(GOEXT)

$(RESULT_DIR)/$1/e_vec_$(2)_$(3)/: $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA)
	mkdir -p $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/

$(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA): $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(GOEXT) instrumentation/benchmark/cpu_vector_benchmark.cpp
	$(GPP) -march=znver2 -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
                            -DDIMENSION=3 \
                            instrumentation/benchmark/cpu_vector_benchmark.cpp \
                            $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT) \
                            -o $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA)

$(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX): instrumentation/benchmark/cpu_vector_benchmark.cpp
	$(GPP) -mcpu=$(XCPU) --target=$(XTRIPLE) -O3 $(CFLAGS) -DLATTICE_$(word 1,$(subst _, ,$(1))) \
                            -DDIMENSION=3 \
                            instrumentation/benchmark/cpu_vector_benchmark.cpp \
                            $(LLVM_DIR)/$(1)/e_vec_$(2)_$(3)/$(LOEXT) \
                            -o $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX)

$(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA).csv: $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA) $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/
	mkdir -p $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/
	$(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA) > $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA).csv

$(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX).csv: $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX) $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/
	mkdir -p $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/
	$(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX) > $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX).csv

CPU_EXP += $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA)
CPU_EXP_RES += $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TA).csv
XCPU_EXP += $(LLVM_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX)
XCPU_EXP_RES += $(RESULT_DIR)/$1/e_vec_$(2)_$(3)/benchmark_$(TX).csv
endef


##############################################################
# Setup all the experiments
##############################################################
# $(eval $(call SS_GPU_RENDER,cc_zp_element,45,1))
# $(eval $(call SS_GPU_RENDER,fcc_v2,4,4))
# $(eval $(call SS_GPU_RENDER,cc_tpcubic,60,60))
# $(eval $(call SS_GPU_RENDER,bcc_v2,4,1))
# $(eval $(call SS_GPU_RENDER,cc_fcc_v3,28,28))
# $(eval $(call SS_GPU_RENDER,cc_bcc_v3,47,11))
# $(eval $(call SS_GPU_RENDER,bcc_tp_lin,15,6))
# $(eval $(call SS_GPU_RENDER,bcc_rdodecq,3,1))
# $(eval $(call SS_GPU_RENDER,bcc_tp_quad,33,26))
# $(eval $(call SS_GPU_RENDER,cc_6dir,12,1))
# $(eval $(call SS_GPU_RENDER,fcc_tp_quad,43,14))
# $(eval $(call SS_GPU_RENDER,cc_rdodecl,15,8))
# $(eval $(call SS_GPU_RENDER,bcc_v3,16,1))
# $(eval $(call SS_GPU_RENDER,fcc_6dir,5,5))
# $(eval $(call SS_GPU_RENDER,fcc_v3,2,1))
# $(eval $(call SS_GPU_RENDER,fcc_tp_lin,26,26))
# $(eval $(call SS_GPU_RENDER,bcc_rdodecl,1,1))
# $(eval $(call SS_GPU_RENDER,cc_toctq,33,26))
# $(eval $(call SS_GPU_RENDER,bcc_toctq,24,1))
# $(eval $(call SS_GPU_RENDER,cc_tplin,7,6))
# $(eval $(call SS_GPU_RENDER,cc_rdodecq,34,1))
# $(eval $(call SS_GPU_RENDER,cc_bcc_v2,32,31))
# $(eval $(call SS_GPU_RENDER,cc_fcc_v2,8,1))
# $(eval $(call SS_GPU_RENDER,cc_tpquad,11,1))
#
include experiments/*.experiments.mk
include experiments_cc/*.experiments.mk
include experiments_bcc/*.experiments.mk
# include experiments/fcc_6dir_cpu.experiments.mk

# Missing experiments
# include experiments/fcc_v2_cpu.experiments.mk
# include experiments/fcc_tp_lin_cpu.experiments.mk
# include experiments/fcc_tp_quad_cpu.experiments.mk
# include experiments_cc/cc_6dir_cpu.experiments.mk

# $(eval $(call SS_CPU_VEC_EXPERIMENT,fcc_v2,3,2))
# include experiments_cc/cc_tpquad_gpu.experiments.mk
# $(eval $(call SS_GPU_LIN_EXPERIMENT,cc_tplin,4))

build_cpu_val: $(CPU_VAL)
cpu_validate: $(CPU_VAL_RES)

build_xcpu_val: $(XCPU_VAL)
xcpu_validate: $(XCPU_VAL_RES)

build_cpu_expr: $(CPU_EXP)
run_cpu_expr: $(CPU_EXP_RES)

build_xcpu_expr: $(XCPU_EXP)
run_xcpu_expr: $(XCPU_EXP_RES)

build_gpu_expr: $(GPU_EXP)
run_gpu_expr: $(GPU_EXP_RES)

build_gpu_lin_expr: $(GPU_LIN_EXP)
run_gpu_lin_expr: $(GPU_LIN_EXP_RES)

build_volren_expr: $(GPU_REN_EXP)
run_volren_expr: $(GPU_REN_RES)