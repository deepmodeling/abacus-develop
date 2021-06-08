#==========================
# VARS
#==========================
FORTRAN        = gfortran
CPLUSPLUS      = g++
CPLUSPLUS_MPI  = mpicxx

LAPACK_DIR     = /usr/local
FFTW_DIR       = /usr/local
BOOST_DIR      = /usr/local
ELPA_DIR       = /usr/local
CEREAL_DIR     = /usr/local
SCALAPACK_DIR  = /usr/local

OBJ_DIR = obj
NP      = 4

#==========================
# LIB and INCLUDE
#==========================
LAPACK_INCLUDE_DIR = ${LAPACK_DIR}/include
LAPACK_LIB_DIR     = ${LAPACK_DIR}/lib

HONG_FFTW        = -D__FFTW3
FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_LIB_DIR     = ${FFTW_DIR}/lib
FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

ELPA_LIB_DIR = ${ELPA_DIR}/lib
ELPA_INCLUDE_DIR = ${ELPA_DIR}/include/elpa-2018.11.001
ELPA_LIB     = -L${ELPA_LIB_DIR} -lelpa -Wl,-rpath=${ELPA_LIB_DIR}

LIBXC_INCLUDE_DIR = ${LIBXC_DIR}/include
LIBXC_LIB_DIR     = ${LIBXC_DIR}/lib
LIBXC_LIB         = -L${LIBXC_LIB_DIR} -lxc -Wl,-rpath=${LIBXC_LIB_DIR}

CEREAL_INCLUDE_DIR = ${CEREAL_DIR}/include
#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lgfortran -lm -openmp -lpthread ${SCALAPACK_DIR}/lib/libscalapack.a ${LAPACK_LIB_DIR}/libopenblas.a ${FFTW_LIB} ${ELPA_LIB}
INCLUDES = -I. -Icommands -I${LAPACK_INCLUDE_DIR} -I${FFTW_INCLUDE_DIR} -I${LIBXC_INCLUDE_DIR} -I${CEREAL_INCLUDE_DIR} -I${ELPA_INCLUDE_DIR}
#==========================
# OPTIMIZE OPTIONS
#==========================
OPTS     = ${INCLUDES} -O2 -std=c++11 -march=native -Wall -fopenmp
