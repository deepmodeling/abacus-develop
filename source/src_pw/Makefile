# This is the Makefile of ABACUS-PW API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = mpiicpc
CUDA_COMPILE = nvcc
OBJ_DIR = pw_obj
NP      = 12

##==========================
## FFTW package needed 
##==========================
#FFTW_DIR = /home/qianrui/intelcompile/impi_fftw
#FFTW_DIR = /home/mohan/1_Software/impi_fftw-3.3.8
FFTW_DIR = /home/haosen/intel-source/fftw-3.3.4
HONG_FFTW        = -D__FFTW3
#FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_INCLUDE_DIR = ${MKLROOT}/include/fftw
#FFTW_LIB_DIR     = ${FFTW_DIR}/lib
#FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

##==========================
## CUDA needed 
##==========================
CUDA_DIR = /usr/local/cuda-11.0
CUDA_INCLUDE_DIR	= ${CUDA_DIR}/include 
CUDA_LIB_DIR		= ${CUDA_DIR}/lib64
CUDA_LIB 			= -L${CUDA_LIB_DIR} -lcufft -lcublas -lcudart


#==========================
# LAPACK package needed 
#==========================
LAPACK_DIR    		= $(MKLROOT)
LAPACK_INCLUDE_DIR	= ${LAPACK_DIR}/include 
LAPACK_LIB_DIR		= ${LAPACK_DIR}/lib/intel64
LAPACK_LIB 			= -L${LAPACK_LIB_DIR} -lmkl_intel_lp64 -lmkl_intel_thread\
					  -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 #-Wl,--end-group -Wl,-rpath=${LAPACK_LIB_DIR}

##==========================
## LIBXC package needed 
##==========================
#LIBXC_DIR    		= /home/qianrui/intelcompile/intel_libxc4.3.4
#LIBXC_INCLUDE_DIR	= ${LIBXC_DIR}/include 
#LIBXC_LIB_DIR		= ${LIBXC_DIR}/lib
#LIBXC_LIB 			= -L${LIBXC_LIB_DIR} -Wl,-rpath=${LIBXC_LIB_DIR} -lxc

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB} ${LAPACK_LIB} ${CUDA_LIB}
CUDA_LIBS = ${CUDA_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR} -I${LAPACK_INCLUDE_DIR} -I${CUDA_INCLUDE_DIR}
CUDA_INCLUDES = -I. -Icommands -I${CUDA_INCLUDE_DIR}

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g
OPTS_CUDA = ${INCLUDES} -std=c++11 

include Makefile.Objects

VPATH=../src_global\
:../module_base\
:../module_cell\
:../module_symmetry\
:../src_parallel\
:../src_io\
:../module_relaxation\
:../module_md\
:../module_symmetry\
:../\

#==========================
# Define HONG
#==========================
#HONG= -DMETIS -DMKL_ILP64 -D__NORMAL ${HONG_FFTW} 
#HONG= ${HONG_FFTW}  -D__MPI -DUSE_LIBXC
HONG= ${HONG_FFTW}  -D__MPI -D__NOMD -D__CUDA
FP_OBJS_0=main.o\
$(OBJS_BASE)\
$(OBJS_PW)\
$(OBJS_PARALLEL)\
${OBJS_IO}\
${OBJS_CELL}\
$(OBJS_MAIN)\
${OBJS_IONS}

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
pw : 
	@ make init
	@ make -j $(NP) parallel

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

parallel : ${FP_OBJS} 
	${CPLUSPLUS_MPI} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/main.o:../main.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} ../main.cpp -o ${OBJ_DIR}/main.o
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

# rules of .cu files
${OBJ_DIR}/%.o:%.cu
	${CUDA_COMPILE} ${OPTS_CUDA} ${LIBS} -I/home/haosen/intel/compilers_and_libraries_2017.8.262/linux/mpi/intel64/include -c ${HONG} $< -o $@


.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}.x ]; then rm -f ${VERSION}.x; fi
