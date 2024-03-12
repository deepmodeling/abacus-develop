# This is the Makefile of ABACUS-ORB API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = mpiicpc 
FFTW_DIR = /home/wenfei/codes/FFTW
OBJ_DIR = obj_deepks
NP      =  4

#==========================
# FFTW package needed 
#==========================
HONG_FFTW        = -D__FFTW3 -D__DEEPKS
FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_LIB_DIR     = ${FFTW_DIR}/lib
FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

#==========================
# libtorch and libnpy
#==========================
LIBTORCH_DIR = /home/wenfei/codes/libtorch
LIBNPY_DIR = /home/wenfei/codes/libnpy

LIBTORCH_INCLUDE_DIR =  -isystem ${LIBTORCH_DIR}/include -isystem ${LIBTORCH_DIR}/include/torch/csrc/api/include
LIBTORCH_LIB_DIR= ${LIBTORCH_DIR}/lib
LIBTORCH_LIB = -L${LIBTORCH_LIB_DIR} -ltorch -lc10 -Wl,-rpath,${LIBTORCH_LIB_DIR} -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch_cpu.so" -Wl,--as-needed ${LIBTORCH_LIB_DIR}/libc10.so -lpthread -Wl,--no-as-needed,"${LIBTORCH_LIB_DIR}/libtorch.so" -Wl,--as-needed 

CNPY_INCLUDE_DIR = ${LIBNPY_DIR}

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB} ${LIBTORCH_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR} ${LIBTORCH_INCLUDE_DIR} -I${CNPY_INCLUDE_DIR} 

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++14 -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../../module_base\
:../../module_io\
:../../module_hamilt_pw/hamilt_pwdft\
:../../module_basis/module_ao\
:../../module_neighbor\
:../../module_cell\
:../../\
:../\
:./\

#==========================
# Define HONG
#==========================
HONG= -DMETIS -DMKL_ILP64 -D__LCAO ${HONG_FFTW}

FP_OBJS_0=main.o\
LCAO_deepks_test.o\
LCAO_deepks_test_prep.o\
$(OBJS_MAIN)\
$(OBJS_IO)\
$(OBJS_BASE)\
$(OBJS_CELL)\
$(OBJS_ORB)\
$(OBJS_NEIGHBOR)\
$(OBJS_PW)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
DEEPKS : 
	@ make init
	@ make -j $(NP) serial

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

serial : ${FP_OBJS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
