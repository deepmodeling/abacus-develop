CXX=mpicxx \
    cmake -B build \
        -DUSE_DSP=ON \
        -DENABLE_LCAO=OFF \
        -DFFTW3_DIR=/vol8/appsoftware/fftw/ \
        -DFFTW3_LIBRARY=/vol8/appsoftware/fftw/lib/libfftw3.so \
        -DFFTW3_OMP_LIBRARY=/vol8/appsoftware/fftw/lib/libfftw3_omp.so \
        -DFFTW3_FLOAT_LIBRARY=/vol8/appsoftware/fftw/lib/libfftw3f.so \
        -DLAPACK_DIR=/vol8/appsoftware/openblas/0.3.21/lib \
        -DDIR_MTBLAS_LIBRARY=/vol8/home/dptech_zyz1/develop/packages/libmtblas_abacus.so