#include "basic_funcs.h"
#include "spin_constrain.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::print_termination()
{
    print_2d("after-optimization spin: (print in the inner loop): ", this->Mi_);
    std::cout << "Inner optimization for lambda ends." << std::endl;
    std::cout << "===============================================================================" << std::endl;
}

template <>
bool SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::check_rms_stop(int outer_step, int i_step, double rms_error)
{
    std::cout << "Step (Outer -- Inner) =  " << outer_step << " -- " << std::setw(5) << i_step + 1
              << "       RMS =" << rms_error << std::endl;
    if (rms_error < this->sc_thr_ || i_step == this->nsc_ - 1)
    {
        if (rms_error < this->sc_thr_)
        {
            std::cout << "Meet convergence criterion ( < " << this->sc_thr_ << " ), exit." << std::endl;
        }
        else if (i_step == this->nsc_ - 1)
        {
            std::cout << "Reach maximum number of steps ( " << this->nsc_ << " ), exit." << std::endl;
        }
        this->print_termination();
        return true;
    }
    return false;
}