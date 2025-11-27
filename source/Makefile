# This is the Makefile of ABACUS
ABACUS_ROOT := $(dir $(firstword $(MAKEFILE_LIST)))

BUILD_DIR ?= build
OBJ_DIR = $(BUILD_DIR)/obj
BIN_DIR = $(BUILD_DIR)/bin
BUILD_INFO_DIR = $(BUILD_DIR)/source_io

include $(ABACUS_ROOT)Makefile.vars

#==========================
# Compiler information 
#==========================
INCLUDES = -I$(ABACUS_ROOT) -I$(ABACUS_ROOT)source_main -I$(ABACUS_ROOT)source_main/commands -I$(ABACUS_ROOT)commands -I$(BUILD_INFO_DIR) -I$(ABACUS_ROOT)../ -I$(ABACUS_ROOT)source_base/module_container

LIBS = -lm -lpthread
OPTS = -std=c++14 -pedantic -m64 ${INCLUDES}
HONG = -D__LCAO
HONG += -D__ELPA

# An FFT-based spherical Bessel transform algorithm has become
# the default since v3.3.4. Comment out the following line to use
# the original Simpson-based spherical Bessel transform.
HONG += -DUSE_NEW_TWO_CENTER

ifeq ($(OPENMP), ON)
    ELPA_NAME = elpa_openmp
else
    ELPA_NAME = elpa
endif

SOURCE_DIR = $(ABACUS_ROOT)
ifeq ($(findstring mpi, $(CXX)), mpi)
     # We do not support EXX in sequential version temporarily.
     HONG += -D__MPI  -DUSE_CEREAL_SERIALIZATION
     MPI = ON
     TESTNP = 4
     suffix = mpi
     ifeq ($(findstring mpii, $(CXX)), mpii)
         INTEL=ON
     endif
else
     TESTNP = 1
     suffix = serial
     ifeq ($(findstring i, $(CXX)), i)
         INTEL=ON
     endif
endif

ifeq ($(DEBUG), ON)
    INTEL = OFF
    HONG += -DDEBUG
    ifeq ($(MPI), ON)
        CXX = mpicxx
    else
        CXX = g++
    endif
    OPTS += -O0 -fsanitize=address -fno-omit-frame-pointer -Wall -g #It can check segmental defaults
else
    # FIXME -Ofast is not compatible with the new two-center integration algorithm
    HONG += -O3 -march=native -DNDEBUG
endif

ifeq ($(INTEL), ON)
    OPTS += -lifcore  -qopenmp
    ##==========================
    ## MKL package
    ##==========================
    HONG  += -D__FFTW3
    LIBS += -L${MKLROOT}/lib/intel64 -Wl,--start-group -lmkl_intel_lp64\
	       -lmkl_intel_thread -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -Wl,-rpath=${MKLROOT}/lib/intel64
    INCLUDES += -I${MKLROOT}/include/fftw

    #==========================
    # ELPA package
    #==========================
    ELPA_LIB_DIR = ${ELPA_DIR}/lib
    LIBS     += -L${ELPA_LIB_DIR} -l${ELPA_NAME}  -Wl,-rpath=${ELPA_LIB_DIR}
    INCLUDES += -I${ELPA_INCLUDE_DIR}
else
    OPTS += -fpermissive
    ifeq ($(OPENMP), ON)
        OPTS += -fopenmp 
    endif
    ##==========================
    ## FFTW package
    ##==========================
    FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
    FFTW_LIB_DIR     = ${FFTW_DIR}/lib
    HONG  += -D__FFTW3
    LIBS += -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}
    INCLUDES += -I${FFTW_INCLUDE_DIR}
    
    #==========================
    # OPENBLAS, SCALAPACK, ELPA
    #==========================
    LIBS     += -L${OPENBLAS_LIB_DIR} -L${SCALAPACK_LIB_DIR} -L${ELPA_DIR}/lib -l${ELPA_NAME} -lscalapack -lopenblas -lgfortran\
                -Wl,-rpath=${ELPA_DIR}/lib -Wl,-rpath=${SCALAPACK_LIB_DIR} -Wl,-rpath=${OPENBLAS_LIB_DIR}
    INCLUDES += -I${ELPA_INCLUDE_DIR}
endif


#==========================
# Cereal package
#==========================
INCLUDES += -I${CEREAL_DIR}/include



##==========================
## CUDA needed 
##==========================
ifeq ($(GPU), CUDA)
    CUDA_COMPILE = nvcc
    CUDA_DIR            = /usr/local/cuda-11.0
    CUDA_INCLUDE_DIR	= ${CUDA_DIR}/include 
    CUDA_LIB_DIR		= ${CUDA_DIR}/lib64
    CUDA_LIB 			= -L${CUDA_LIB_DIR} -lcufft -lcublas -lcudart
    HONG += -D__CUDA
    OPTS_CUDA = ${INCLUDES} -std=c++11 
    INCLUDES += -I${CUDA_INCLUDE_DIR}
    LIBS += ${CUDA_LIB}
endif

ifdef LIBXC_DIR
    ##==========================
    ## LIBXC package 
    ##==========================
    LIBXC_INCLUDE_DIR	= ${LIBXC_DIR}/include 
    LIBXC_LIB_DIR		= ${LIBXC_DIR}/lib
    HONG                += -DUSE_LIBXC
    INCLUDES            += -I${LIBXC_INCLUDE_DIR}
    LIBS 			    += -L${LIBXC_LIB_DIR} -Wl,-rpath=${LIBXC_LIB_DIR} -lxc
    PWTAG               += 'LIBXC_DIR=${LIBXC_DIR}'
endif

ifdef LIBRI_DIR
    INCLUDES += -I${LIBRI_DIR}/include
    INCLUDES += -I${LIBCOMM_DIR}/include
    HONG                += -DUSE_LIBRI
    HONG += -D__EXX -DEXX_H_COMM=2 -DEXX_DM=3 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1
endif

ifdef LIBNPY_DIR
    CNPY_INCLUDE_DIR = -I${LIBNPY_DIR}/include
    ifeq ($(ENABLE_MLALGO), ON)
        HONG += -D__MLALGO
        # HONG += -D__MLALGO
    endif
    INCLUDES += $(CNPY_INCLUDE_DIR)
endif

ifdef LIBTORCH_DIR
    LIBTORCH_INCLUDE_DIR =  -I${LIBTORCH_DIR}/include -I${LIBTORCH_DIR}/include/torch/csrc/api/include
    INCLUDES += $(LIBTORCH_INCLUDE_DIR)
    LIBTORCH_LIB_DIR= ${LIBTORCH_DIR}/lib
    LIBS += -L${LIBTORCH_LIB_DIR} -ltorch -lc10 -Wl,-rpath,${LIBTORCH_LIB_DIR} -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch_cpu.so" -Wl,--as-needed ${LIBTORCH_LIB_DIR}/libc10.so -lpthread -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch.so" -Wl,--as-needed 
endif

ifdef DeePMD_DIR
    HONG  += -D__DPMD -DHIGH_PREC 
    OPTS  += -Wl,--no-as-needed

    ##==========================
    ## DeePMD-kit package
    ##==========================
    DeeP_LIB_DIR = ${DeePMD_DIR}/lib
    DeeP_INCLUDE_DIR = ${DeePMD_DIR}/include
    INCLUDES += -I${DeeP_INCLUDE_DIR}
    ifneq ($(wildcard ${DeePMD_DIR}/include/deepmd/deepmd.hpp), )
        HONG  += -D__DPMDC
        LIBS += -L${DeeP_LIB_DIR} -ldeepmd_c -Wl,-rpath=${DeeP_LIB_DIR}
    else
        LIBS += -L${DeeP_LIB_DIR} -ldeepmd_cc -Wl,-rpath=${DeeP_LIB_DIR}
    endif
endif

ifdef TensorFlow_DIR
    #==========================
    # TensorFlow package
    #==========================
    TensorFlow_LIB_DIR = ${TensorFlow_DIR}/lib
    LIBS += -L${TensorFlow_LIB_DIR} -ltensorflow_cc -Wl,-rpath=${TensorFlow_LIB_DIR}
    TensorFlow_INCLUDE_DIR = ${TensorFlow_DIR}/include
    INCLUDES += -I${TensorFlow_INCLUDE_DIR}
endif

ifdef PEXSI_DIR
    OBJS_ABACUS += ${OBJS_HSOLVER_PEXSI}
    INCLUDES += -I${PEXSI_DIR}/include -I${PARMETIS_DIR}/include -I${DSUPERLU_DIR}/include
    LIBS += -L${PEXSI_DIR}/lib -lpexsi -L${DSUPERLU_DIR}/lib -lsuperlu_dist -L${PARMETIS_DIR}/lib -lparmetis -lmetis
    HONG += -D__PEXSI
endif

include $(ABACUS_ROOT)Makefile.Objects
VPATH := $(patsubst ./%,$(ABACUS_ROOT)%,$(VPATH))
HEADERS := $(patsubst %,$(ABACUS_ROOT)%,$(HEADERS))

#==========================
# Optional HONG
#==========================
HONG += -D__SELINV -DMETIS 
# ifeq ($(MEM_CHECK), ON)
#     HONG += -D_MCD_CHECK -DWIN32 -DMCD_VERBOSE
# endif

#==========================
# OBJECTS NEEDED
#==========================

ifeq ($(MPI), ON)
    ifdef LIBRI_DIR
        OBJS_ABACUS += $(OBJS_MODULE_RI) ${OBJS_SRC_RI}
    endif
endif

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OBJS_ABACUS})

###### START of ABACUS INFO PART ######

# Rule to generate the build_info.h file using the shell script.
# This rule is executed whenever build_info.h.in, the script, or any of the
# relevant Makefile variables change.
$(BUILD_INFO_DIR)/build_info.h:$(ABACUS_ROOT)source_io/build_info.h.in $(ABACUS_ROOT)../generate_build_info.sh
	@echo "Generating build information header..."
	@# Make sure directory is available
	@mkdir -p $(dir $@)
	@# Export all necessary Makefile variables as environment for the script
	@export BUILD_TYPE="Release"; \
	if [ "$(DEBUG)" = "ON" ]; then export BUILD_TYPE="Debug"; export ENABLE_ASAN="yes"; else export ENABLE_ASAN="no"; fi; \
	export CXX="${CXX}"; \
	export CXXFLAGS="${OPTS} ${HONG}"; \
	export MKLROOT="${MKLROOT}"; \
	export ELPA_DIR="${ELPA_DIR}"; \
	export FFTW_DIR="${FFTW_DIR}"; \
	export CEREAL_DIR="${CEREAL_DIR}"; \
	export LIBXC_DIR="${LIBXC_DIR}"; \
	export LIBRI_DIR="${LIBRI_DIR}"; \
	export DeePMD_DIR="${DeePMD_DIR}"; \
	export TensorFlow_DIR="${TensorFlow_DIR}"; \
	export PEXSI_DIR="${PEXSI_DIR}"; \
	export ENABLE_MPI="no"; \
	if [ "$(MPI)" = "ON" ]; then export ENABLE_MPI="ON"; fi; \
	export USE_CUDA="no"; \
	if [ "$(GPU)" = "CUDA" ]; then export USE_CUDA="ON"; fi; \
	$(ABACUS_ROOT)../generate_build_info.sh $@
	@echo "Build information header generated."

    # Specific rule for parse_args.o to enforce dependency on build_info.h.
    # This overrides the generic %.o:%.cpp rule for this specific file.
${OBJ_DIR}/parse_args.o: $(ABACUS_ROOT)source_io/parse_args.cpp $(BUILD_INFO_DIR)/build_info.h
	@mkdir -p $(dir $@)
	${CXX} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

###### END of ABACUS INFO PART ######

#==========================
# MAKING OPTIONS
#==========================
abacus:
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	@ $(MAKE) $(BIN_DIR)/${VERSION}.$(suffix)

test:
	@ $(MAKE) abacus
	@ cd $(ABACUS_ROOT)../tests/integrate/;sh Autotest.sh -a $(realpath $(BIN_DIR))/ABACUS.mpi -n $(TESTNP)

pw $(BIN_DIR)/${VERSION}-PW.x:
	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	@ cd $(ABACUS_ROOT)source_pw/module_pwdft; $(MAKE) CXX=${CXX} GPU=${GPU} DEBUG=$(DEBUG) FFTW_DIR=$(FFTW_DIR) OPENBLAS_LIB_DIR=$(OPENBLAS_LIB_DIR) ${PWTAG}
	@ cp $(ABACUS_ROOT)source_pw/module_pwdft/${VERSION}-PW.x $(realpath $(BIN_DIR))/${VERSION}-PW.x

$(BIN_DIR)/${VERSION}.$(suffix) : ${FP_OBJS} ${PDIAG_OBJS} ${HEADERS}
	@mkdir -p $(dir $@)
	${CXX} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${LIBS} -o  $(BIN_DIR)/${VERSION}.$(suffix)

#==========================
# rules
#==========================
# Note: The specific rule for parse_args.o above is more precise.
# This generic rule will apply to all other .cpp files.
${OBJ_DIR}/%.o:$(realpath $(BIN_DIR))%.cpp
	${CXX} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
	
.PHONY:clean test
clean:
	@ echo "Cleaning build directory: $(BUILD_DIR)"
	@ rm -rf $(BUILD_DIR)
