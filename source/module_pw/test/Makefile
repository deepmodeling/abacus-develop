#======================================================================
# Please set
# e.g. make CC=mpiicpc or make CC=icpc
#======================================================================
CC = mpiicpc
# mpiicpc:   compile intel parallel version
# icpc:      compile intel serial version
# mpicxx:    compile gnu parallel version
# g++:       compile gnu serial version
GPU = OFF
# OFF:  do not use GPU
# CUDA: use CUDA 
DEBUG = OFF
# ON:   use gnu compiler and check segmental defaults
# OFF:  nothing
FLOAT = OFF
# ON: compile both float and double lib
# OFF: compile double lib
GTEST_DIR = /home/qianrui/gnucompile/g_gtest
# Must give gtest directory
FFTW_DIR = /home/qianrui/gnucompile/fftw_3.3.8
# directory of fftw package. only used when CC=mpicxx/g++
#======================================================================


#==========================
# Compiler information 
#==========================
HONG = -D__NORMAL
INCLUDES = -I.
LIBS = 
OPTS = -Ofast -march=native -std=c++11 -m64  ${INCLUDES}
OBJ_DIR = obj
LIBNAME = libpw.a
GEN     = OFF
TIME    = OFF
ifeq ($(findstring mpi, $(CC)), mpi)
     HONG += -D__MPI
     MPI = ON
     ifeq ($(findstring mpii, $(CC)), mpii)
         INTEL=ON
     endif
else
     ifeq ($(findstring i, $(CC)), i)
         INTEL=ON
     endif
endif


ifeq ($(DEBUG), ON)
    INTEL = OFF
    ifeq ($(MPI), ON)
        CC = mpicxx
    else
        CC = g++
    endif
    HONG += -D__DEBUG
    OPTS += -fsanitize=address -fno-omit-frame-pointer -Wall -g #It can check segmental defaults
endif

ifeq ($(INTEL), ON)
##==========================
## MKL
##==========================
    LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
    INCLUDES += -I${MKLROOT}/include/fftw
    OPTS += -simd -qopenmp
else
##==========================
## fftw 
##==========================
    FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
    FFTW_LIB_DIR     = ${FFTW_DIR}/lib
    LIBS += -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}
    INCLUDES += -I${FFTW_INCLUDE_DIR}
    OPTS += -fopenmp
endif

ifeq ($(FLOAT), ON)
    HONG += -D__MIX_PRECISION
    ifneq ($(INTEL), ON)
        LIBS += -lfftw3f
    endif
endif

##==========================
## CUDA 
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
## GTEST 
##==========================
GTESTOPTS = -I${GTEST_DIR}/include -L${GTEST_DIR}/lib -lgtest -lpthread



#==========================
# Objects
#==========================
VPATH=../../src_parallel\
:../../module_base\
:../

MATH_OBJS0=matrix.o\
matrix3.o\
tool_quit.o\
mymath3.o\
timer.o\
pw_basis.o\
pw_distributer.o\
pw_init.o\
pw_transform.o\
pw_distributeg.o\
pw_distributeg_method1.o\
pw_distributeg_method2.o\
fft.o\
pw_basis_k.o\
pw_transform_k.o

OTHER_OBJS0=global_variable.o\
parallel_global.o

TESTFILE0=test1-1-1.o\
test1-1-2.o\
test1-2.o\
test1-2-2.o\
test1-3.o\
test1-4.o\
test1-5.o\
test2-1-1.o\
test2-1-2.o\
test2-2.o\
test2-3.o\
test3-1.o\
test3-2.o\
test3-3.o\
test4-1.o\
test4-2.o\
test4-3.o\
test4-4.o\
test4-5.o

MATH_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${MATH_OBJS0})
OTHER_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OTHER_OBJS0})
TESTFILE=$(patsubst %.o, ${OBJ_DIR}/%.o, ${TESTFILE0})


#==========================
# MAKING OPTIONS
#==========================
pw_test.exe: ${LIBNAME} ${OTHER_OBJS} ${TESTFILE} pw_test.cpp test_tool.cpp
	${CC} ${OPTS} ${HONG} pw_test.cpp test_tool.cpp ${TESTFILE} ${LIBNAME} ${OTHER_OBJS} ${LIBS} -o pw_test.exe ${GTESTOPTS}
	@ if [ $(GEN) == "ON" ]; then $(MAKE) gen.exe ; fi
	@ if [ $(TIME) == "ON" ]; then $(MAKE) time.exe ; fi

lib $(LIBNAME): $(MATH_OBJS)
	ar -rcv $(LIBNAME) $(MATH_OBJS)

gen.exe: ${LIBNAME} ${OTHER_OBJS} generate.cpp test_tool.cpp
	${CC} ${OPTS} ${HONG} generate.cpp test_tool.cpp $(LIBNAME) ${OTHER_OBJS} ${LIBS} -o gen.exe

time.exe: ${LIBNAME} time.cpp
	${CC} ${OPTS} ${HONG} time.cpp $(LIBNAME)  ${LIBS} -o time.exe

${OBJ_DIR}/%.o:%.cpp $(OBJ_DIR)/readme.log
	${CC} ${OPTS} -c ${HONG} $< -o $@ ${GTESTOPTS}

$(OBJ_DIR)/readme.log:
	@if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@echo "This is a temporary dirctory to store obj files." > $(OBJ_DIR)/readme.log


.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e $(LIBNAME) ]; then rm -f $(LIBNAME); fi
	@ if [ -e pw_test.exe ]; then rm -f pw_test.exe; fi
	@ if [ -e gen.exe ]; then rm -f gen.exe; fi
	@ if [ -e time.exe ]; then rm -f time.exe; fi
