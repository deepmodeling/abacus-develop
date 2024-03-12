# This is the Makefile of ABACUS-PW API
#======================================================================
# Users set
# e.g. make INTEL=ON or make INTEL=OFF GPU=OFF
#======================================================================
CXX = mpiicpc
# mpiicpc:   compile intel parallel version
# mpicxx:    compile gnu parallel version
# icpc:      compile gnu sequential version
# g++:       compile gnu sequential version
# We do not support sequential version with g++ or icpc yet.
GPU = OFF    
# OFF:  do not use GPU
# CUDA: use CUDA 
DEBUG = OFF
# ON:   use gnu compiler and check segmental defaults
# OFF:  nothing
FFTW_DIR = /home/qianrui/gnucompile/fftw_3.3.8
# directory of fftw package. only used when INTEL = OFF
OPENBLAS_LIB_DIR   = /home/qianrui/gnucompile/openblas/lib
# directory of librefblas.a and liblapack.a. Only used when INTEL = OFF
#======================================================================


#==========================
# Compiler information 
#==========================
INCLUDES = -I. -Icommands -I../../
LIBS = -lm -lpthread
OPTS = -Ofast -march=native -std=c++11 -m64 ${INCLUDES}
ifeq ($(findstring mpi, $(CXX)), mpi)
     HONG += -D__MPI -D__EXX -DEXX_H_COMM=2 -DUSE_CEREAL_SERIALIZATION -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1
     MPI = ON
     ifeq ($(findstring mpii, $(CXX)), mpii)
         INTEL=ON
     endif
else
     ifeq ($(findstring i, $(CXX)), i)
         INTEL=ON
     endif
endif

ifeq ($(DEBUG), ON)
    INTEL = OFF
    HONG += -D__DEBUG
    ifeq ($(MPI), ON)
        CXX = mpicxx
    else
        CXX = g++
    endif
    OPTS +=  -fsanitize=address -fno-omit-frame-pointer -Wall -g #It can check segmental defaults
endif
ifeq ($(INTEL), ON)
    OPTS += -qopenmp -lifcore
    ##==========================
    ## MKL package
    ##==========================
    HONG  += -D__FFTW3
    LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread\
    				  -lmkl_core
    INCLUDES += -I${MKLROOT}/include/fftw
    
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
    # OPENBLAS package
    #==========================
    LIBS += -L${OPENBLAS_LIB_DIR}  -lopenblas -lgfortran
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



include ../../Makefile.Objects
VPATH:=$(addprefix ../../, ${VPATH})
OBJ_DIR = pw_obj

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OBJS_ABACUS_PW})

#==========================
# MAKING OPTIONS
#==========================
pw : 
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	$(MAKE) ${VERSION}-PW.x

${VERSION}-PW.x : ${FP_OBJS} 
	${CXX} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}-PW.x 
#==========================
# rules
#==========================
${OBJ_DIR}/main.o:../main.cpp
	${CXX} ${OPTS} -c ${HONG} ../main.cpp -o ${OBJ_DIR}/main.o
	
${OBJ_DIR}/%.o:%.cpp
	${CXX} ${OPTS} -c ${HONG} $< -o $@

# rules of .cu files
${OBJ_DIR}/%.o:%.cu
	${CUDA_COMPILE} ${OPTS_CUDA} ${LIBS} -I/home/haosen/intel/compilers_and_libraries_2017.8.262/linux/mpi/intel64/include -c ${HONG} $< -o $@

	

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}-PW.x ]; then rm -f ${VERSION}-PW.x; fi