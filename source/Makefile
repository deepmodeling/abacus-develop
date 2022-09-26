# This is the Makefile of ABACUS
include Makefile.vars
#==========================
# Compiler information 
#==========================
INCLUDES = -I. -Icommands -I../
LIBS = -lifcore -lm -lpthread
OPTS = -Ofast -march=native -std=c++14 -pedantic -m64 ${INCLUDES}
HONG = -D__LCAO
OBJ_DIR = obj
BIN_DIR = ../bin
ifeq ($(findstring mpi, $(CC)), mpi)
     # We do not support EXX in serial version temporarily.
     HONG += -D__MPI -D__EXX -DEXX_H_COMM=2 -DUSE_CEREAL_SERIALIZATION -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1
     MPI = ON
     TESTNP = 4
     suffix = mpi
     ifeq ($(findstring mpii, $(CC)), mpii)
         INTEL=ON
     endif
else
     TESTNP = 1
     suffix = serial
     ifeq ($(findstring i, $(CC)), i)
         INTEL=ON
     endif
endif

ifeq ($(DEBUG), ON)
    INTEL = OFF
    HONG += -D__DEBUG
    ifeq ($(MPI), ON)
        CC = mpicxx
    else
        CC = g++
    endif
    OPTS +=  -fsanitize=address -fno-omit-frame-pointer -Wall -g #It can check segmental defaults
endif

ifeq ($(INTEL), ON)
    OPTS += -simd -qopenmp
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
    LIBS     += -L${ELPA_LIB_DIR} -lelpa  -Wl,-rpath=${ELPA_LIB_DIR}
    INCLUDES += -I${ELPA_INCLUDE_DIR}
else
    OPTS += -fopenmp -fpermissive
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
    LIBS     += -L${OPENBLAS_LIB_DIR} -L${SCALAPACK_LIB_DIR} -L${ELPA_DIR}/lib -lelpa -lscalapack -lopenblas -lgfortran\
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
FP_OBJS_0=$(OBJS_FIRST_PRINCIPLES)\
main.o\

ifeq ($(MPI), ON)
    FP_OBJS_0 += $(OBJS_PDIAG) $(PDIAG_MR_0) $(OBJS_EXX)
endif

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

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
	@ cd src_pw; $(MAKE) CC=${CC} GPU=${GPU} DEBUG=$(DEBUG) FFTW_DIR=$(FFTW_DIR) OPENBLAS_LIB_DIR=$(OPENBLAS_LIB_DIR)
	@ cp src_pw/${VERSION}-PW.x $(BIN_DIR)/${VERSION}-PW.x

$(BIN_DIR)/${VERSION}.$(suffix) : ${FP_OBJS} ${PDIAG_OBJS} ${PDIAG_MR} ${HEADERS}
	${CC} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${PDIAG_MR} ${LIBS} -o  $(BIN_DIR)/${VERSION}.$(suffix)

# ${VERSION}.fp_mem.x : ${FP_OBJS} ../src_parallel/mcd.o ${HEADERS}
# 	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
# 	${CC} $(FP_OBJS) ../src_parallel/mcd.o ${LIBS} -o ${VERSION}.fp_mem.x
# 	@ cp ${VERSION}.fp_mem.x $(BIN_DIR)/${VERSION}.fp_mem.x

#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CC} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
	
.PHONY:clean test
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -d $(BIN_DIR) ]; then rm -rf $(BIN_DIR); fi
	@ if [ -e ${VERSION}.* ]; then rm -f ${VERSION}.*; fi
	@ cd src_pw; make clean
# @ if [ -e ${VERSION}.fp_mem.x ]; then rm -f ${VERSION}.fp_mem.x; fi
