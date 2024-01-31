SERVER=FUGAKU # for FUGAKU
#SERVER=LOCAL # for gfortarn
#SERVER=LOCAL_ifx # for ifx

# -Dlittle using little endian
# -Dart_ro using momentum for the artifical viscosity
# -Dentropy using entrop

PPC := $(PPC) -Dposixio_write
PPC := $(PPC) -Dposixio_read
PPC := $(PPC) -Dremap_2d_assign

########################################################
### gfortran
ifeq ($(SERVER), LOCAL)
FC = /usr/bin/mpif90
FFLAGS0 = -cpp -Dlittle $(PPC) -mcmodel=medium -fallow-argument-mismatch
FFLAGS = $(FFLAGS0) -O3
endif

#######################################################################
ifeq ($(SERVER), LOCAL_ifx)
FC = mpiifx
FFLAGS0 = -fpp -Dlittle -Difort -shared-intel -mcmodel=large $(PPC)
FFLAGS = $(FFLAGS0) -O3
endif

#######################################################################
### FUGAKU
ifeq ($(SERVER), FUGAKU)
FC = mpifrtpx # FUGAKU
FFLAGS := -Cpp -Dlittle -Dfujitsu $(PPC) -Kfast,openmp,parallel,array_private,ocl -Koptmsg=2 -Nlst=t -fPIC
FFLAGS := $(FFLAGS) -x- -xstmt_no,100
endif

# Source directory
# Exec directory
EXE_DIR = ..

UPDATE = update

SRC_FILES := $(wildcard ../src/all/*.F90) 

OBJ_FILES := $(notdir $(SRC_FILES:.F90=.o))
SRC_DIR := $(dir $(SRC_FILES))
VPATH := $(SRC_DIR)

$(UPDATE) : $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $(EXE_DIR)/a.out $(OBJ_FILES) $(FFLAGS_AFTER) && \
	touch $(UPDATE)

%.o: %.F90
	$(FC) $(FFLAGS) -c $< $(FFLAGS_AFTER)

clean:
	rm -rf $(MAIN) $(OBJ_FILES) $(EXE_DIR)/a.out *.mod *.lst *.i $(UPDATE)

## This dependnecy files is created with F90_deps.py in main directory
include R2D2.deps
