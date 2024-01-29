# This is the Makefile of ABACUS
include Makefile.vars
#==========================
# Compiler information 
#==========================

INCLUDES = -I. -Icommands -I../ -Imodule_base/module_container

LIBS = -lm -lpthread
OPTS = -std=c++14 -pedantic -m64 ${INCLUDES}
HONG = -D__LCAO
HONG += -D__ELPA
ifeq ($(OPENMP), ON)
    ELPA_NAME = elpa_openmp
else
    ELPA_NAME = elpa
endif
OBJ_DIR = obj
BIN_DIR = ../bin
SOURCE_DIR = .
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
    HONG += -Ofast -march=native -DNDEBUG
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

ifdef LIBTORCH_DIR
 ifdef LIBNPY_DIR
      LIBTORCH_INCLUDE_DIR =  -I${LIBTORCH_DIR}/include -I${LIBTORCH_DIR}/include/torch/csrc/api/include
      LIBTORCH_LIB_DIR= ${LIBTORCH_DIR}/lib
      LIBTORCH_LIB = -L${LIBTORCH_LIB_DIR} -ltorch -lc10 -Wl,-rpath,${LIBTORCH_LIB_DIR} -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch_cpu.so" -Wl,--as-needed ${LIBTORCH_LIB_DIR}/libc10.so -lpthread -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch.so" -Wl,--as-needed 
      CNPY_INCLUDE_DIR = -I${LIBNPY_DIR}/include
      HONG += -D__DEEPKS
      INCLUDES += $(LIBTORCH_INCLUDE_DIR) $(CNPY_INCLUDE_DIR)
      LIBS += $(LIBTORCH_LIB)
  endif
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

    #==========================
    # TensorFlow package
    #==========================
    ifndef TensorFlow_DIR
        TensorFlow_DIR = ${DeePMD_DIR}
    endif
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

include Makefile.Objects

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

#==========================
# MAKING OPTIONS
#==========================
abacus:
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	@ $(MAKE) $(BIN_DIR)/${VERSION}.$(suffix)

test:
	@ $(MAKE) abacus
	@ cd ../tests/integrate/;sh Autotest.sh -a ../../../bin/ABACUS.mpi -n $(TESTNP)

pw $(BIN_DIR)/${VERSION}-PW.x:
	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	@ cd module_hamilt_pw/hamilt_pwdft; $(MAKE) CXX=${CXX} GPU=${GPU} DEBUG=$(DEBUG) FFTW_DIR=$(FFTW_DIR) OPENBLAS_LIB_DIR=$(OPENBLAS_LIB_DIR) ${PWTAG}
	@ cp module_hamilt_pw/hamilt_pwdft/${VERSION}-PW.x $(BIN_DIR)/${VERSION}-PW.x

$(BIN_DIR)/${VERSION}.$(suffix) : ${FP_OBJS} ${PDIAG_OBJS} ${HEADERS}
	${CXX} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${LIBS} -o  $(BIN_DIR)/${VERSION}.$(suffix)


#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CXX} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
	
.PHONY:clean test
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -d $(BIN_DIR) ]; then rm -rf $(BIN_DIR); fi
	@ cd module_hamilt_pw/hamilt_pwdft; make clean
