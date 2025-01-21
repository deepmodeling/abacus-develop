#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "module_base/scalapack_connector.h" // Cpxgemr2d
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_psi/psi.h"

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

template <typename Device = base_device::DEVICE_CPU>
class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO<std::complex<double>, double>
{
  public:
    ESolver_KS_LCAO_TDDFT();

    ~ESolver_KS_LCAO_TDDFT();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

  protected:
    virtual void hamilt2density_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    virtual void update_pot(UnitCell& ucell, const int istep, const int iter) override;

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter) override;

    virtual void after_scf(UnitCell& ucell, const int istep) override;

    //! wave functions of last time step
    psi::Psi<std::complex<double>>* psi_laststep = nullptr;

    //! Hamiltonian of last time step
    std::complex<double>** Hk_laststep = nullptr;

    //! Overlap matrix of last time step
    std::complex<double>** Sk_laststep = nullptr;

    const int td_htype = 1;

    //! Control heterogeneous computing of the TDDFT solver
    bool use_tensor = false;
    bool use_lapack = false;

  private:
    void weight_dm_rho();
};

} // namespace ModuleESolver
#endif
