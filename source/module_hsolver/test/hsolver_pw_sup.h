#include "module_basis/module_pw/pw_basis_k.h"

namespace ModulePW
{

PW_Basis::PW_Basis(){};
PW_Basis::~PW_Basis(){};

void PW_Basis::initgrids(
    const double lat0_in, //unit length (unit in bohr)
    const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors (unit in lat0)
    const double gridecut //unit in Ry, ecut to set up grids
)
{
    return;
}

void PW_Basis::initgrids(
    const double lat0_in,
    const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
    const int nx_in, int ny_in, int nz_in
)
{
    return;
}

void PW_Basis::distribute_r()
{
    return;
}

PW_Basis_K::PW_Basis_K()
{
    this->nks = 1;
    this->npwk_max = 3;
    //used for update_precondition
    this->gk2 = new double[3];
    this->npwk = new int[1];
    this->npwk[0] = 3;
    this->tpiba2 = 1.0;
}

PW_Basis_K::~PW_Basis_K()
{
    delete[] this->gk2;
    delete[] this->npwk;
}

double& PW_Basis_K::getgk2(const int ik, const int igl) const
{
    this->gk2[igl] = (ik + igl) * 1.5;
    return this->gk2[igl];
}

FFT::FFT()
{
}

FFT::~FFT()
{
}

}//namespace ModulePW

#include "module_hsolver/diago_cg.h"
#include "module_hsolver/diago_david.h"
#include "module_hsolver/diago_iter_assist.h"

template<typename T> consts<T>::consts() {}
template class consts<std::complex<float>>;
template class consts<std::complex<double>>;

namespace hsolver
{

template<typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(
    const std::string& basis_type,
    const std::string& calculation)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template<typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(
    const std::string& basis_type,
    const std::string& calculation,
    const bool& need_subspace,
    const Func& subspace_func,
    const Real& pw_diag_thr,
    const int& pw_diag_nmax,
    const int& nproc_in_pool)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    need_subspace_ = need_subspace;
    subspace_func_ = subspace_func;
    pw_diag_thr_ = pw_diag_thr;
    pw_diag_nmax_ = pw_diag_nmax;
    nproc_in_pool_ = nproc_in_pool;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template<typename T, typename Device>
DiagoCG<T, Device>::~DiagoCG()
{
    delete this->one_;
    delete this->zero_;
    delete this->neg_one_;
}

template<typename T, typename Device>
void DiagoCG<T, Device>::diag(
    const Func& hpsi_func, 
    const Func& spsi_func, 
    ct::Tensor& psi,
    ct::Tensor& eigen,
    const ct::Tensor& prec)
{
    auto n_bands = psi.shape().dim_size(0);
    auto n_basis = psi.shape().dim_size(1);
    auto psi_pack = psi.accessor<T, 2>();
    auto eigen_pack = eigen.accessor<Real, 1>();
    //do something
    for(int ib = 0;ib<n_bands;ib++)
    {
        eigen_pack[ib] = 0.0;
        for(int ig = 0;ig<n_basis;ig++)
        {
            psi_pack[ib][ig] += T(2.0, 0.0);
            eigen_pack[ib] += psi_pack[ib][ig].real();
        }
        eigen_pack[ib] /= n_basis;
    }
    DiagoIterAssist<T, Device>::avg_iter += 1.0;
    return;
}

template class DiagoCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoCG<std::complex<double>, base_device::DEVICE_CPU>;

template <typename T, typename Device>
DiagoDavid<T, Device>::DiagoDavid(const Real* precondition_in, 
                                  int diago_david_ndim_in,
                                  bool use_paw_in,
                                  const diag_comm_info& diag_comm_in)
    : david_ndim(diago_david_ndim_in), use_paw(use_paw_in), diag_comm(diag_comm_in)
{
    this->device = base_device::get_device_type<Device>(this->ctx);
    this->precondition = precondition_in;

    this->one = &this->cs.one;
    this->zero = &this->cs.zero;
    this->neg_one = &this->cs.neg_one;

    test_david = 2;
    // 1: check which function is called and which step is executed
    // 2: check the eigenvalues of the result of each iteration
    // 3: check the eigenvalues and errors of the last result
    // default: no check
}

template <typename T, typename Device>
DiagoDavid<T, Device>::~DiagoDavid()
{
    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->sphi);
    delmem_complex_op()(this->ctx, this->hcc);
    delmem_complex_op()(this->ctx, this->scc);
    delmem_complex_op()(this->ctx, this->vcc);
    delmem_complex_op()(this->ctx, this->lagrange_matrix);
    base_device::memory::delete_memory_op<Real, base_device::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue);
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::diag(hamilt::Hamilt<T, Device>* phm_in,
                                      psi::Psi<T, Device>& psi,
                                      Real* eigenvalue_in)
{
    //do something
    for(int ib = 0;ib<psi.get_nbands();ib++)
    {
        eigenvalue_in[ib] = 0.0;
        for(int ig = 0;ig<psi.get_nbasis();ig++)
        {
            psi(ib, ig) += T(1.0, 0.0);
            eigenvalue_in[ib] += psi(ib, ig).real();
        }
        eigenvalue_in[ib] /= psi.get_nbasis();
    }
    DiagoIterAssist<T, Device>::avg_iter += 1.0;
    return;
}
template class DiagoDavid<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoDavid<std::complex<double>, base_device::DEVICE_CPU>;

template class DiagoIterAssist<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>;

}//namespace hsolver

#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
namespace hamilt
{

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_CPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<float>, base_device::DEVICE_CPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>* phm_in)
{
    for (int i = 0; i < wvf.size(); i++)
    {
        wvf.get_pointer()[i] = std::complex<float>((float)i + 1, 0);
    }
}

template <>
void diago_PAO_in_pw_k2(const base_device::DEVICE_CPU* ctx,
                        const int& ik,
                        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* phm_in)
{
    for (int i = 0; i < wvf.size(); i++)
    {
        wvf.get_pointer()[i] = std::complex<double>( (double)i+1, 0) ;
    }
}

}//namespace hsolver

template class hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU>;
template class hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU>;