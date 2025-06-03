//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================

#include "driver.h"
#include "fftw3.h"
#include "module_base/parallel_global.h"
#include "module_io/parse_args.h"
#include "module_parameter/parameter.h"
#include "module_basis/module_pw/module_fft/fft_bundle.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "module_basis/module_pw/module_fft/test.cuh"
int main(int argc, char** argv)
{
    /*
    read the arguement in the command-line,
    with "abacus -v", the program exit and returns version info,
    with no arguments, the program continues.
    */
       std::cout << "FFT Bundle Example" << std::endl;

    // Example usage of make_unique
    ModulePW::test();
    // Note: The following lines are commented out as they require specific FFT implementations
    // Uncomment and implement the FFT operations as needed
    
    // auto fft_bundle = make_unique<FFT_Bundle>("gpu", "single");

    // Initialize FFT parameters
    // fft_bundle->initfft(256, 256, 256, 64, 64, 1, 1, 1, false, false, false);

    // Perform FFT operations
    // std::complex<float> input[256];
    // std::complex<float> output[256];
    // fft_bundle->fftxyfor(input, output);

    std::cout << "FFT operation completed." << std::endl;

    return 0;
}