#======================================================================
# Please set
# e.g. make CXX=mpiicpc or make CXX=icpc
#======================================================================
CXX = mpiicpc
# mpiicpc:   compile intel parallel version
# icpc:      compile intel sequential version
# mpicxx:    compile gnu parallel version
# g++:       compile gnu sequential version
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
#======================================================================


#==========================
# Compiler information 
#==========================
HONG = -D__NORMAL
INCLUDES = -I. -I../../../ -I../../../source_main -I../../../source_base/module_container
LIBS = 
OPTS = -Ofast -march=native -std=c++11 -m64  ${INCLUDES}
OBJ_DIR = obj
LIBNAME = libpw.a
GEN     = OFF
TIME    = OFF
ifeq ($(findstring mpi, $(CXX)), mpi)
     HONG += -D__MPI
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
    ifeq ($(MPI), ON)
        CXX = mpicxx
    else
        CXX = g++
    endif
    HONG += -D__DEBUG
    OPTS += -fsanitize=address -fno-omit-frame-pointer -static-libasan -Wall -g #It can check segmental defaults
endif

ifeq ($(INTEL), ON)
##==========================
## MKL
##==========================
    LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
    INCLUDES += -I${MKLROOT}/include/fftw
    OPTS +=  -qopenmp
else
##==========================
## fftw 
##==========================
    LIBS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64   -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm
    INCLUDES += -I${MKLROOT}/include/fftw
    OPTS += 
endif

ifeq ($(FLOAT), ON)
    HONG += -D__ENABLE_FLOAT_FFTW
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
GTESTOPTS = -I${GTEST_DIR}/include -L${GTEST_DIR}/lib -lgtest -lpthread -w



#==========================
# Objects
#==========================
VPATH=../../../source_main\
../../../source_base\
../../../source_base/module_device\
../../../source_base/module_container/ATen/core\
../../../source_base/module_container/ATen\
../../../source_io/module_parameter\
../module_fft\
../\

MATH_OBJS0=matrix.o\
matrix3.o\
tool_quit.o\
mymath.o\
timer.o\
pw_basis.o\
pw_distributer.o\
pw_init.o\
pw_transform.o\
pw_distributeg.o\
pw_distg_method1.o\
pw_distg_method2.o\
pw_basis_k.o\
pw_basis_sup.o\
pw_transform_k.o\
memory.o\
memory_op.o\
depend_mock.o\
parameter.o\
fft_cpu.o\
fft_cpu_float.o\
fft_bundle.o\



OTHER_OBJS0=


TESTFILE0=test1_1_1.o\
test1_1_2.o\
test1_2.o\
test1_2_2.o\
test1_3.o\
test1_4.o\
test1_5.o\
test2_1_1.o\
test2_1_2.o\
test2_2.o\
test2_3.o\
test3_1.o\
test3_2.o\
test3_3.o\
test3_3_2.o\
test4_1.o\
test4_2.o\
test4_3.o\
test4_4.o\
test4_5.o\
test5_1_1.o\
test5_1_2.o\
test5_2_1.o\
test5_2_2.o\
test5_3_1.o\
test5_4_1.o\
test5_4_2.o\
test6_1_1.o\
test6_1_2.o\
test6_2_1.o\
test6_2_2.o\
test6_3_1.o\
test6_4_1.o\
test6_4_2.o\
test7_1.o\
test7_2_1.o\
test7_3_1.o\
test7_3_2.o\
test8_1.o\
test8_2_1.o\
test8_3_1.o\
test8_3_2.o\
test_big.o\
test_other.o

MATH_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${MATH_OBJS0})
OTHER_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OTHER_OBJS0})
TESTFILE=$(patsubst %.o, ${OBJ_DIR}/%.o, ${TESTFILE0})


#==========================
# MAKING OPTIONS
#==========================
pw_test.exe: ${LIBNAME} ${OTHER_OBJS} ${TESTFILE} pw_test.cpp test_tool.cpp
	${CXX} ${OPTS} ${HONG} pw_test.cpp test_tool.cpp ${TESTFILE} ${LIBNAME} ${OTHER_OBJS} ${LIBS} -o pw_test.exe ${GTESTOPTS}
	@ if [ $(GEN) == "ON" ]; then $(MAKE) gen.exe ; fi
	@ if [ $(TIME) == "ON" ]; then $(MAKE) time.exe ; fi

lib $(LIBNAME): $(MATH_OBJS)
	ar -rcv $(LIBNAME) $(MATH_OBJS)

gen.exe: ${LIBNAME} ${OTHER_OBJS} generate.cpp test_tool.cpp
	${CXX} ${OPTS} ${HONG} generate.cpp test_tool.cpp $(LIBNAME) ${OTHER_OBJS} ${LIBS} -o gen.exe

time.exe: ${LIBNAME} time.cpp
	${CXX} ${OPTS} ${HONG} time.cpp $(LIBNAME)  ${LIBS} -o time.exe

${OBJ_DIR}/%.o:%.cpp $(OBJ_DIR)/readme.log
	${CXX} ${OPTS} -c ${HONG} $< -o $@ ${GTESTOPTS}

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
