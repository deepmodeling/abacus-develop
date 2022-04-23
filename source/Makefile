# This is the Makefile of ABACUS

ifdef GNU
    include Makefile.gnu
else
    include Makefile.system
endif
include Makefile.Objects

VPATH=./src_global\
:./module_orbital\
:./module_neighbor\
:./module_symmetry\
:./module_cell\
:./module_base\
:./module_md\
:./module_surchem\
:./module_deepks\
:./module_xc\
:./module_esolver\
:./src_pw\
:./src_lcao\
:./src_ions\
:./src_io\
:./src_parallel\
:./src_pdiag\
:./src_pdiag/MRRR\
:./src_ri\
:./\


#==========================
# Define HONG
#==========================
#(1)no mpi
HONG_SER = -D__FP -D__LCAO  ${HONG_FFTW}
HONG_SER_SELINV = -D__FP -D__LCAO ${HONG_FFTW} -D__SELINV
HONG_GDB = -g -D__FP -D__LCAO  ${HONG_FFTW}
#(2)mpi
HONG_MPI = -D__FP -D__LCAO  ${HONG_FFTW} -D__MPI
# mohan comment out 2021-02-06, add -DUSE_LIBXC=0 if you want to use LIBXC
#HONG_MPI_SELINV = -D__FP ${HONG_FFTW} ${HONG_LAPACK} -D__MPI -D__SELINV -DMETIS -DEXX_DM=3 -DEXX_H_COMM=2 -DUSE_LIBXC=0 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1 -DUSE_CEREAL_SERIALIZATION
HONG_MPI_SELINV = -D__FP ${HONG_FFTW} ${HONG_LAPACK} -D__MPI -D__LCAO -D__SELINV -DMETIS -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1 -DUSE_CEREAL_SERIALIZATION
HONG_MPI_SELINV_20210523 = -D__FP ${HONG_FFTW} ${HONG_LAPACK} -D__LCAO -D__MPI -D__SELINV -DMETIS -DEXX_DM=3 -DEXX_H_COMM=2 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1 -DUSE_CEREAL_SERIALIZATION -D__EXX
HONG_DEEPKS = ${HONG_MPI_SELINV_20210523} -D__DEEPKS	#caoyu add 2021-07-15 , use it in DeePKS. Add LIBTORCH_DIR and LIBNPY_DIR, modify ${LIBS} and build with std=c++14 !!
#(3)memory
#(3)memory
HONG_MEM = ${HONG_FFTW} -D__FP -D_MCD_CHECK -DWIN32 -DMCD_VERBOSE

#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!
#! change series version or parallel version~~~
HONG=${HONG_MPI_SELINV_20210523}
#HONG=${HONG_SER_SELINV}
#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!

#==========================
# OPTIMIZE OPTIONS
#==========================
#(1)gdb
#xiaohui modify 2013-12-27
OPTS_GDB = -g -W -I${BOOST_INCLUDE_DIR}
#(2)gprof
#POTS = -pg

#==========================
# OBJECTS NEEDED
#==========================
FP_OBJS_0=$(OBJS_FIRST_PRINCIPLES)\
$(OBJS_PARALLEL)\
$(OBJS_TOOLS) \
$(OBJS_COMMON)\
$(OBJS_SURCHEM)\
main.o\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})
FP_OBJS_serial=$(patsubst %.o, ${OBJ_DIR_serial}/%.o, ${FP_OBJS_0})
PDIAG_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OBJS_PDIAG})
PDIAG_MR=$(patsubst %.o, ${OBJ_DIR}/%.o, ${PDIAG_MR_0})

#==========================
# MAKING OPTIONS
#==========================
fp_mpi : 
	@ make init
	@ make -j $(NP) parallel

fp : 
	@ make init_serial
	@ make -j $(NP) serial

fp_gbd : 
	@ make init
	@ make -j $(NP) serial_gbd

fp_gprof : 
	@ make init
	@ make -j $(NP) serial_gprof

fp_mem : 
	@ make init
	@ make -j $(NP) serial_mem

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi
	@ if [ ! -d ../bin ]; then mkdir ../bin; fi

init_serial :
	@ if [ ! -d $(OBJ_DIR_serial) ]; then mkdir $(OBJ_DIR_serial); fi
	@ if [ ! -d $(OBJ_DIR_serial)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR_serial)/README; fi
	@ if [ ! -d ../bin ]; then mkdir ../bin; fi

parallel : ${FP_OBJS} ${PDIAG_OBJS} ${PDIAG_MR} ${HEADERS}
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${PDIAG_MR} ${LIBS} -o  ../bin/${VERSION}.mpi

serial : ${FP_OBJS_serial} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS_serial) ${LIBS} -o ${VERSION}.fp.x 

serial_gdb : ${FP_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS_GDB} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp_gdb.x

serial_gprof : ${FP_OBJS} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp_gprof.x 

serial_mem : ${FP_OBJS} ../src_parallel/mcd.o ${HEADERS}
	${CPLUSPLUS} $(FP_OBJS) ../src_parallel/mcd.o ${LIBS} -o ${VERSION}.fp_mem.x

#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
${OBJ_DIR}/%.o:%.f
	${FORTRAN} -c ${HONG} $< -o $@	 

${OBJ_DIR_serial}/%.o:%.cpp
	${CPLUSPLUS} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
${OBJ_DIR_serial}/%.o:%.f
	${FORTRAN} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -d $(OBJ_DIR_serial) ]; then rm -rf $(OBJ_DIR_serial); fi
	@ if [ -d ../bin ]; then rm -rf ../bin; fi
	@ cd src_pw; make clean
