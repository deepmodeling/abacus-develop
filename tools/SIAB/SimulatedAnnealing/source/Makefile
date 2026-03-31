#===========================================================
# AUTHOR : mohan
# DATE START FROM : 2009-04-14
# LAST UPDATE : 2025-08-28
# FUNCTION : spillage makefile
#===========================================================
MKLROOT=/home/mohan/intel/oneapi/mkl/021.3.0/

#------------------------------
# for serial version: make s
#------------------------------
CPLUSPLUS = icpc 

#------------------------------
# for parallel version: make p
#------------------------------
#CPLUSPLUS = mpiicpc

CXXFLAGS = -std=c++11 -O3

#CXXFLAGS = -std=c++11 -O2 -xCORE-AVX2 -align -fma -ftz -fomit-frame-pointer -finline-functions \
           -I${MKLROOT}/include -I${MKLROOT}/include/intel64/ilp64 \
           -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lpthread -limf  

LIBS=-L${MKLROOT}/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm

TOOLS_OBJS=src_tools/mathzone.o \
			src_tools/matrix.o \
			src_tools/matrix3.o \
			src_tools/realarray.o \
			src_tools/intarray.o \
			src_tools/complexmatrix.o \
			src_tools/complexarray.o\
			src_tools/timer.o \
			src_tools/inverse_matrix.o\
			src_tools/inverse_matrix_iterate.o\
			src_tools/Simulated_Annealing.o\

SPILLAGE_OBJS=src_spillage/main.o \
			  src_spillage/read_INPUT.o\
			  src_spillage/ReadData.o\
			  src_spillage/Calculate_C4.o\
			  src_spillage/Coefficients.o\
			  src_spillage/MultiZeta.o\
			  src_spillage/SpillageStep.o\
			  src_spillage/SpillageValue.o\
			  src_spillage/Step_Data.o\
			  src_spillage/Inverse_Matrix_S.o\
			  src_spillage/Orthogonal.o\
			  src_spillage/Type_Information.o\
			  src_spillage/Metropolis.o\
			  src_spillage/Metropolis_2.o\
			  src_spillage/Out_Orbital.o\
			  src_spillage/Plot_Psi.o\
			  src_spillage/Psi_Second.o\
			  src_spillage/tools.o \
			  src_spillage/common.o \

PARALLEL_OBJS=src_parallel/parallel_global.o\
	src_parallel/parallel_common.o\
	src_parallel/parallel_reduce.o\
	src_parallel/parallel_kpoints.o	

PW_OBJS=src_pw/pw_basis.o\
	src_pw/pw_complement.o\
	src_pw/heapsort.o\
	src_pw/numerical_basis.o\
	src_pw/bessel_basis.o\
	src_pw/ylm_real.o\
	src_pw/memory_calculation.o\
			  
OBJS = ${TOOLS_OBJS} \
	   ${SPILLAGE_OBJS} \
	   ${PARALLEL_OBJS} \
	   ${PW_OBJS}\

help : 
	@echo 'make target list: '
	@echo 's : Series version of Monte Carlo Spillage Program.'
	@echo 'p : Parallel Version of Monte Carlo Spillage Program.'
	@echo 'note : some notes.'

s : ${OBJS} 
	${CPLUSPLUS} ${CXXFLAGS} -o SIA_s.exe ${OBJS} ${LIBS}

p : ${OBJS}
	${CPLUSPLUS} -D__MPI ${CXXFLAGS} -o SIA_p.exe ${OBJS} ${LIBS}

note :
	@echo 'atom type label           : it '
	@echo 'atom label for each type  : ia '
	@echo 'The angular momentum      : l '
	@echo 'The multiplicity          : n '
	@echo 'The magnetic quantum num  : m '
	@echo 'The k point index         : ik '
	@echo 'The band index            : ib '
	@echo 'The eigenvalue of Jlq     : ie '
	@echo 'The radius cutoff (a.u.)  : rcut '
	@echo 'The state is              : psi(ib,ik) '
	@echo 'The local orbitals is     : phi(it,ia,l,n,m) '
	@echo 'The index of local basis  : nwfc(it,l,n)'
	@echo 'The index of local basis  : nwfc2(it,ia,l,n,m) '
	@echo 'The Spherical Bessel func : Jl(ie,rcut) '
	@echo 'The SB transform          : \int Jl(ie,rcut) * Jl'
	@echo 'The index of Jlq coef     : Coef(it,l,n,ie)'
	@echo 'The local orbitals formu  : phi=\sum_{ie} Coef(it,l,n,ie) * Jl(ie,rcut)'
	@echo 'The spillage formula is   : S = <psi|(1-P)|psi>'
	@echo 'The projector P is        : P = \sum_{mu,nu} |phi_mu> S^{-1} <phi_nu| '
	@echo 'The orthogonal formula is : |psi2> = (1-P1) |psi1>'
	@echo 'So we can fonud that      : P1|psi2> = 0 '
	@echo 'the Q matrix is           : Q = <Jlq|psi> '
	@echo 'the S matrix is           : S = <Jlq|Jlq> '
	@echo 'start index of band       : BANDS_START(count from 1, not 0) '
	@echo 'end index of band         : BANDS_END '

.cpp.o:
	${CPLUSPLUS} ${CXXFLAGS} -c $< -o $@

.PHONY:clean
clean:
	${RM} *.o *.exe core* src_tools/*.o src_spillage/*.o tests_s/ORBITAL* tests_s/ORBITAL* src_parallel/*.o src_pw/*.o
