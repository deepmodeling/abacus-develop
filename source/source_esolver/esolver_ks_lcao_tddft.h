#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "source_base/module_container/ATen/core/tensor.h"   // ct::Tensor
#include "source_base/module_external/scalapack_connector.h" // Cpxgemr2d
#include "source_lcao/module_rt/td_info.h"
#include "source_lcao/module_rt/velocity_op.h"

namespace ModuleESolver
{
//------------------------ MPI gathering and distributing functions ------------------------//
// This struct is used for collecting matrices from all processes to root process
template <typename T>
struct Matrix_g
{
    std::shared_ptr<T> p;
    size_t row;
    size_t col;
    std::shared_ptr<int> desc;
};

// Collect matrices from all processes to root process
template <typename T>
void gatherMatrix(const int myid, const int root_proc, const hamilt::MatrixBlock<T>& mat_l, Matrix_g<T>& mat_g)
{
    const int* desca = mat_l.desc; // Obtain the descriptor of the local matrix
    int ctxt = desca[1];           // BLACS context
    int nrows = desca[2];          // Global matrix row number
    int ncols = desca[3];          // Global matrix column number

    if (myid == root_proc)
    {
        mat_g.p.reset(new T[nrows * ncols]); // No need to delete[] since it is a shared_ptr
    }
    else
    {
        mat_g.p.reset(new T[nrows * ncols]); // Placeholder for non-root processes
    }

    // Set the descriptor of the global matrix
    mat_g.desc.reset(new int[9]{1, ctxt, nrows, ncols, nrows, ncols, 0, 0, nrows});
    mat_g.row = nrows;
    mat_g.col = ncols;

    // Call the Cpxgemr2d function in ScaLAPACK to collect the matrix data
    Cpxgemr2d(nrows, ncols, mat_l.p, 1, 1, const_cast<int*>(desca), mat_g.p.get(), 1, 1, mat_g.desc.get(), ctxt);
}
//------------------------ MPI gathering and distributing functions ------------------------//

template <typename TR, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO<std::complex<double>, TR>
{
  public:
    ESolver_KS_LCAO_TDDFT();

    ~ESolver_KS_LCAO_TDDFT();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

  protected:
    virtual void runner(UnitCell& cell, const int istep) override;

    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    void store_h_s_psi(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver);

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    void print_step();

    //! Wave function for all k-points of last time step
    psi::Psi<std::complex<double>>* psi_laststep = nullptr;

    //! Hamiltonian for all k-points of last time step
    ct::Tensor Hk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Overlap matrix for all k-points of last time step
    ct::Tensor Sk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE);

    //! Control heterogeneous computing of the TDDFT solver
    bool use_tensor = false;
    bool use_lapack = false;

    // Control the device type for Hk_laststep and Sk_laststep
    // Set to CPU temporarily, should wait for further GPU development
    static constexpr ct::DeviceType ct_device_type_hs = ct::DeviceType::CpuDevice;

    //! Total steps for evolving the wave function
    int totstep = -1;

    //! Velocity matrix for calculating current
    Velocity_op<TR>* velocity_mat = nullptr;

    TD_info* td_p = nullptr;

    //! Restart flag
    bool restart_done = false;

  private:
    void weight_dm_rho(const UnitCell& ucell);
};

} // namespace ModuleESolver
#endif // ESOLVER_KS_LCAO_TDDFT_H
