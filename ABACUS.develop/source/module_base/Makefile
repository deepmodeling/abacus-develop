# This is the Makefile of ABACUS-base API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = icpc 
FFTW_DIR = /home/mohan/1_Software/impi_fftw-3.3.8
OBJ_DIR = obj
NP      = 1

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands 

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../src_global\
:./\

#==========================
# Define HONG
#==========================
HONG= -DMETIS -DMKL_ILP64 -D__NORMAL

FP_OBJS_0=main.o\
$(OBJS_BASE)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
base : 
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
