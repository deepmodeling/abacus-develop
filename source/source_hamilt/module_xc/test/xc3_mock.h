// I'm mocking FFT here because it is not possible to write
// unit tests with FFT
#include <fstream>

namespace ModulePW
{
namespace
{
// Preserve the shared prefix and use zero for mock outputs without a
// one-to-one counterpart between real- and reciprocal-space buffers.
template <typename FPTYPE, typename InputType>
void mock_real2recip(const InputType* in, std::complex<FPTYPE>* out, const int nrxx, const int nrecip, const bool add, const FPTYPE factor)
{
    for (int i = 0; i < nrecip; ++i)
    {
        const std::complex<FPTYPE> value = (i < nrxx) ? std::complex<FPTYPE>(in[i]) : std::complex<FPTYPE>(0.0, 0.0);
        if (add)
        {
            out[i] += factor * value;
        }
        else
        {
            out[i] = value;
        }
    }
}

template <typename FPTYPE>
void mock_recip2real(const std::complex<FPTYPE>* in,
                     std::complex<FPTYPE>* out,
                     const int nrecip,
                     const int nrxx,
                     const bool add,
                     const FPTYPE factor)
{
    for (int i = 0; i < nrxx; ++i)
    {
        const std::complex<FPTYPE> value = (i < nrecip) ? -ModuleBase::IMAG_UNIT * in[i] : std::complex<FPTYPE>(0.0, 0.0);
        if (add)
        {
            out[i] += factor * value;
        }
        else
        {
            out[i] = value;
        }
    }
}

template <typename FPTYPE>
void mock_recip2real(const std::complex<FPTYPE>* in, FPTYPE* out, const int nrecip, const int nrxx, const bool add, const FPTYPE factor)
{
    for (int i = 0; i < nrxx; ++i)
    {
        const FPTYPE value = (i < nrecip) ? (-ModuleBase::IMAG_UNIT * in[i]).real() : FPTYPE(0.0);
        if (add)
        {
            out[i] += factor * value;
        }
        else
        {
            out[i] = value;
        }
    }
}
} // namespace

PW_Basis::PW_Basis() {};
PW_Basis::~PW_Basis() {};

template <typename FPTYPE>
void PW_Basis::real2recip(const FPTYPE* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
{
    mock_real2recip(in, out, nrxx, npw, add, factor);
}
template void PW_Basis::real2recip<double>(const double* in, std::complex<double>* out, bool add, double factor) const;

template <typename FPTYPE>
void PW_Basis::real2recip(const std::complex<FPTYPE>* in, std::complex<FPTYPE>* out, const bool add, const FPTYPE factor) const
{
    mock_real2recip(in, out, nrxx, npw, add, factor);
}
template void PW_Basis::real2recip<double>(const std::complex<double>* in,
                                           std::complex<double>* out,
                                           const bool add,
                                           const double factor) const;

template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out,
                          const bool add,
                          const FPTYPE factor) const // in:(nz, ns)  ; out(nplane,nx*ny)
{
    mock_recip2real(in, out, npw, nrxx, add, factor);
}
template void PW_Basis::recip2real(const std::complex<double>* in, std::complex<double>* out, const bool add, const double factor) const;

template <typename FPTYPE>
void PW_Basis::recip2real(const std::complex<FPTYPE>* in, FPTYPE* out, const bool add, const FPTYPE factor) const
{
    mock_recip2real(in, out, npw, nrxx, add, factor);
}
template void PW_Basis::recip2real(const std::complex<double>* in, double* out, const bool add, const double factor) const;

template <typename FPTYPE>
void PW_Basis_K::recip2real(const std::complex<FPTYPE>* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const // in:(nz, ns)  ; out(nplane,nx*ny)
{
    mock_recip2real(in, out, npwk[ik], nrxx, add, factor);
}
template void PW_Basis_K::recip2real(const std::complex<double>* in,
                                     std::complex<double>* out,
                                     const int ik,
                                     const bool add,
                                     const double factor) const;

ModuleBase::Vector3<double> PW_Basis_K::getgpluskcar(int, int) const
{
    ModuleBase::Vector3<double> x = {1, 2, 3};
    return x;
}

template <typename FPTYPE, typename Device>
void PW_Basis_K::real_to_recip(const Device* ctx,
                               const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const // in:(nplane,nx*ny)  ; out(nz, ns)
{
    mock_real2recip(in, out, nrxx, npwk[ik], add, factor);
}
template <typename FPTYPE, typename Device>
void PW_Basis_K::recip_to_real(const Device* ctx,
                               const std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out,
                               const int ik,
                               const bool add,
                               const FPTYPE factor) const
{
    mock_recip2real(in, out, npwk[ik], nrxx, add, factor);
}

template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_CPU>(const base_device::DEVICE_CPU* ctx,
                                                                         const std::complex<double>* in,
                                                                         std::complex<double>* out,
                                                                         const int ik,
                                                                         const bool add,
                                                                         const double factor) const;
template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_CPU>(const base_device::DEVICE_CPU* ctx,
                                                                         const std::complex<double>* in,
                                                                         std::complex<double>* out,
                                                                         const int ik,
                                                                         const bool add,
                                                                         const double factor) const;
#if __CUDA || __ROCM
template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_GPU>(const base_device::DEVICE_GPU* ctx,
                                                                         const std::complex<double>* in,
                                                                         std::complex<double>* out,
                                                                         const int ik,
                                                                         const bool add,
                                                                         const double factor) const;

template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_GPU>(const base_device::DEVICE_GPU* ctx,
                                                                         const std::complex<double>* in,
                                                                         std::complex<double>* out,
                                                                         const int ik,
                                                                         const bool add,
                                                                         const double factor) const;
#endif

void PW_Basis::initgrids(double, ModuleBase::Matrix3, double) {};
void PW_Basis::distribute_r() {};
void PW_Basis::initgrids(double, ModuleBase::Matrix3, int, int, int) {};

PW_Basis_K::PW_Basis_K() {};
PW_Basis_K::~PW_Basis_K() {};
} // namespace ModulePW

namespace ModuleBase
{
void WARNING_QUIT(const std::string& file, const std::string& description)
{
    std::cout << " " << file << "  warning : " << description << std::endl;
    exit(1);
}
void WARNING(const std::string& file, const std::string& description) {};

void Matrix3::Identity() {};

IntArray::IntArray(int, int) {};
IntArray::~IntArray() {};

void TITLE(const std::string& class_function_name, bool disable) {};
void TITLE(const std::string& class_name, const std::string& function_name, bool disable) {};

} // namespace ModuleBase


UnitCell::UnitCell() {};
UnitCell::~UnitCell() {};

Charge::Charge() {};
Charge::~Charge() {};

Magnetism::Magnetism() {};
Magnetism::~Magnetism() {};

SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}

namespace unitcell
{
void cal_ux(UnitCell& ucell, const int nspin)
{
    ucell.magnet.lsign_ = false;

    ucell.magnet.ux_[0] = 0;
    ucell.magnet.ux_[1] = 1;
    ucell.magnet.ux_[2] = 2;

    ucell.magnet.lsign_ = true;
};
} // namespace unitcell

namespace Parallel_Reduce
{
/// reduce in all process
template <typename T>
void reduce_all(T& object) {};
template <typename T>
void reduce_all(T* object, const int n) {};
template <typename T>
void reduce_pool(T& object) {};
template <typename T>
void reduce_pool(T* object, const int n) {};

template <>
void Parallel_Reduce::reduce_pool<double>(double& object)
{
#ifdef __MPI
    double swap = object;
    MPI_Allreduce(&swap, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return;
}
template void reduce_all<double>(double& object);
template void reduce_all<double>(double* object, const int n);
template void reduce_pool<float>(float& object);
template void reduce_pool<float>(float* object, const int n);
template void reduce_pool<double>(double* object, const int n);
} // namespace Parallel_Reduce
