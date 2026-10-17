#include <RI/ri/LRI.h>

#if defined(EXPECT_LIBRI_GPU)
#if !defined(__GPU_RI) || !defined(__DDLA_RI) || !defined(DDLA_USE_CUDA)
#error "LibRI GPU usage requirements did not reach the consumer"
#endif
#else
#if defined(__GPU_RI) || defined(__DDLA_RI) || defined(DDLA_USE_CUDA)
#error "CPU LibRI consumer unexpectedly requires a GPU backend"
#endif
#endif

int main()
{
    RI::LRI<int, int, 3, double> lri;
#ifdef EXPECT_LIBRI_GPU
    if (lri.cal_mode != RI::LRI_Cal_Mode::CPU)
    {
        return 1;
    }
    lri.cal_mode = RI::LRI_Cal_Mode::GPU;
#endif
    // The consumer must receive the BLAS dependency through LibRI::LibRI.
    const double data[] = {3.0, 4.0};
    const double norm = RI::Blas_Interface::nrm2(2, data);
    return norm == 5.0 ? 0 : 1;
}
