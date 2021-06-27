# This is the Makefile of ABACUS-neighor API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = icpc
#BOOST_DIR = /home/cyFortneu/soft/boost_1_74_0
BOOST_DIR = /home/mohan/1_Software/impi_boost-1.70.0
#FFTW_DIR = /home/cyFortneu/soft/fftw-3.3.8
FFTW_DIR = /home/mohan/1_Software/impi_fftw-3.3.8

OBJ_DIR = obj
NP      = 1

#==========================
# FFTW package needed 
#==========================
HONG_FFTW        = -D__FFTW3
FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_LIB_DIR     = ${FFTW_DIR}/lib
FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

#==========================
# boost package needed 
#==========================
BOOST_INCLUDE_DIR = ${BOOST_DIR}/include

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
OPTS_GDB = -g -W -I${BOOST_INCLUDE_DIR}

INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR} -I${BOOST_INCLUDE_DIR}

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../src_global\
:./\

#==========================
# Define HONG
#==========================
HONG= -DMETIS -DMKL_ILP64 -D__NORMAL ${HONG_FFTW} 

FP_OBJS_0=main.o\
$(OBJS_NEIGHBOR)\
$(OBJS_BASE)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
neighbor : 
	@ make init
	@ make -j $(NP) serial

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

serial : ${FP_OBJS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

serial_gdb : ${FP_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS_GDB} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp_gdb.x


#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
