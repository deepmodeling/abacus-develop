# This is the Makefile of ABACUS-CMD API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = mpiicpc
OBJ_DIR = md_obj
NP      = 8

##==========================
## FFTW3 
##==========================
HONG_FFTW        = -D__FFTW3
FFTW_INCLUDE_DIR = ${MKLROOT}/include/fftw
#FFTW_DIR = /usr/local/fftw
#FFTW_DIR = /home/qianrui/gnucompile/g_fftw-3.3.8
#FFTW_DIR = /home/dell/packages/fftw3_gnu
#FFTW_DIR = /home/mohan/1_Software/g_fftw3.3.8_serial
#FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
#FFTW_LIB_DIR     = ${FFTW_DIR}/lib
#FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

##============================
## DeePMD-kit (need by DPMD)
##============================
HONG_DPMD = -D__DPMD -DHIGH_PREC 
DeeP_DIR = ${deepmd_root}
DeeP_INCLUDE_DIR = ${DeeP_DIR}/include
DeeP_LIB_DIR = ${DeeP_DIR}/lib
DeeP_LIB = -L${DeeP_LIB_DIR} -ldeepmd_cc -ldeepmd_op -ldeepmd_op_cuda -Wl,-rpath=${DeeP_LIB_DIR}

##============================
## TensorFlow (need by DPMD)
##============================
TensorFlow_DIR = ${tensorflow_root}
TensorFlow_INCLUDE_DIR = ${TensorFlow_DIR}/include
TensorFlow_LIB_DIR = ${TensorFlow_DIR}/lib
TensorFlow_LIB = -L${TensorFlow_LIB_DIR} -ltensorflow_cc -Wl,-rpath=${TensorFlow_LIB_DIR}

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB}
#LIBS = -lm -lpthread -openmp -lstdc++ ${FFTW_LIB}
#LIBS = -lm -lpthread -openmp -lstdc++ ${FFTW_LIB} ${DeeP_LIB} ${TensorFlow_LIB} -Wl,--no-as-needed

INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR}
#INCLUDES = -g -I. -Icommands -I${FFTW_INCLUDE_DIR}
#INCLUDES = -g -I. -Icommands -I${FFTW_INCLUDE_DIR} -I${DeeP_INCLUDE_DIR} -I${TensorFlow_INCLUDE_DIR}

#==========================
# OPTIMIZE OPTIONS
#==========================
# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -xHost -m64 -qopenmp -Werror -Wall -pedantic -g
#OPTS = ${INCLUDES} -Ofast  -std=c++11 -march=native -Wall -fpermissive -fopenmp

include Makefile.Objects

VPATH=../module_base\
:../module_neighbor\
:../module_cell\
:../src_parallel\
:../src_io\
:./\
:../\

#==========================
# Define HONG
#==========================
HONG= ${HONG_FFTW} -D__MPI -DMETIS -DMKL_ILP64 -D__NORMAL -D__CMD -D__CELL
#HONG= ${HONG_FFTW} -D__MPI -DMETIS -DMKL_ILP64 -D__NORMAL -D__CMD -D__CELL ${HONG_DPMD}

FP_OBJS_0=main.o\
$(OBJS_MAIN)\
$(OBJS_BASE)\
$(OBJS_CELL)\
$(OBJS_MD)\
$(OBJS_IO)\
$(OBJS_PARALELL)\
$(OBJS_NEIGHBOR)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
cmd : 
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
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}.x ]; then rm -f ${VERSION}.x; fi
