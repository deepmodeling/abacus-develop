# This is the Makefile of ABACUS-PW API
#======================================================================
# User set
# e.g. make INTEL=ON or make INTEL=OFF GPU=OFF
#======================================================================
INTEL = ON
# ON:   use intel compilers  e.g. icpc/mpiicpc
# OFF:  use gnu compilers  e.g. g++/mpicxx
GPU = OFF    
# OFF:  do not use GPU
# CUDA: use CUDA 
DEBUG = OFF
# ON:   use gnu compiler and check segmental defaults
# OFF:  nothing
FFTW_DIR = /home/qianrui/gnucompile/fftw_3.3.8 
# directory of fftw package. only used when INTEL = OFF
LAPACK_LIB_DIR   = /home/qianrui/gnucompile/lapack_3.9.0/lib 
# directory of librefblas.a and liblapack.a. Only used when INTEL = OFF
#======================================================================


#==========================
# Compiler information 
#==========================
INCLUDES = -I. -Icommands -I../
LIBS = -lifcore -lm -lpthread
OPTS = -Ofast -std=c++11 -Wall -m64 -g ${INCLUDES}
ifeq ($(DEBUG), ON)
    INTEL = OFF
    HONG += -D__DEBUG
    OPTS +=  -fsanitize=address -fno-omit-frame-pointer #It can check segmental defaults
endif
ifeq ($(INTEL), ON)
    CPLUSPLUS     = icpc
    CPLUSPLUS_MPI = mpiicpc
    # -pedantic turns off more extensions and generates more warnings
    # -xHost generates instructions for the highest instruction set available on the compilation host processor
    OPTS += -simd -qopenmp -Werror
    ##==========================
    ## MKL package
    ##==========================
    HONG  += -D__FFTW3
    LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread\
    				  -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
    INCLUDES += -I${MKLROOT}/include/fftw
    
else
    CPLUSPLUS     = g++
    CPLUSPLUS_MPI = mpicxx
    OPTS += -march=native -fopenmp -fpermissive
    ##==========================
    ## FFTW package
    ##==========================
    FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
    FFTW_LIB_DIR     = ${FFTW_DIR}/lib
    LIBS += -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}
    INCLUDES += -I${FFTW_INCLUDE_DIR}
    
    #==========================
    # LAPACK package
    #==========================
    LIBS += -L${LAPACK_LIB_DIR}  -llapack -lrefblas -lgfortran
endif

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


##==========================
## LIBXC package needed 
##==========================
#LIBXC_DIR    		= /home/qianrui/intelcompile/intel_libxc4.3.4
#LIBXC_INCLUDE_DIR	= ${LIBXC_DIR}/include 
#LIBXC_LIB_DIR		= ${LIBXC_DIR}/lib
#HONG              += -DUSE_LIBXC
#INCLUDES          += -I${LIBXC_INCLUDE_DIR}
#LIBS 			   += -L${LIBXC_LIB_DIR} -Wl,-rpath=${LIBXC_LIB_DIR} -lxc



include Makefile.Objects

OBJ_DIR = pw_obj

#==========================
# Define HONG
#==========================
HONG += -D__MPI

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
parallel : ${FP_OBJS} 
	${CPLUSPLUS_MPI} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/main.o:../main.cpp
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ echo "This directory contains all of the .o files." > $(OBJ_DIR)/readme.log
	${CPLUSPLUS_MPI} ${OPTS} -c ${HONG} ../main.cpp -o ${OBJ_DIR}/main.o
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} -c ${HONG} $< -o $@

# rules of .cu files
${OBJ_DIR}/%.o:%.cu
	${CUDA_COMPILE} ${OPTS_CUDA} ${LIBS} -I/home/haosen/intel/compilers_and_libraries_2017.8.262/linux/mpi/intel64/include -c ${HONG} $< -o $@

	

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}.x ]; then rm -f ${VERSION}.x; fi
