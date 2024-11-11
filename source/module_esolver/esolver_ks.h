#ifndef ESOLVER_KS_H
#define ESOLVER_KS_H
#include "esolver_fp.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hsolver/hsolver.h"
#include "module_io/cal_test.h"
#include "module_psi/psi.h"

#ifdef __MPI
#include <mpi.h>
#else
#include <chrono>
#endif
#include <cstring>
#include <fstream>
namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS : public ESolver_FP
{
  public:
    //! Constructor
    ESolver_KS();

    //! Deconstructor
    virtual ~ESolver_KS();

    virtual void before_all_runners(const Input_para& inp, UnitCell& cell) override;

    virtual void runner(const int istep, UnitCell& cell) override;

  protected:
    //! Something to do before SCF iterations.
    virtual void before_scf(const int istep) {};

    virtual void init_after_vc(const Input_para& inp, UnitCell& cell) override; // liuyu add 2023-03-09

    //! Something to do before hamilt2density function in each iter loop.
    virtual void iter_init(const int istep, const int iter);

    //! Something to do after hamilt2density function in each iter loop.
    virtual void iter_finish(const int istep, int& iter);

    // calculate electron density from a specific Hamiltonian with ethr
    virtual void hamilt2density_single(const int istep, const int iter, const double ethr);

    // calculate electron states from a specific Hamiltonian
    virtual void hamilt2estates(const double ethr) {};

    // calculate electron density from a specific Hamiltonian
    void hamilt2density(const int istep, const int iter, const double ethr);

    //! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
    virtual void after_scf(const int istep) override;

    //! <Temporary> It should be replaced by a function in Hamilt Class
    virtual void update_pot(const int istep, const int iter) {};

    //! Hamiltonian
    hamilt::Hamilt<T, Device>* p_hamilt = nullptr;

    ModulePW::PW_Basis_K* pw_wfc = nullptr;

    Charge_Mixing* p_chgmix = nullptr;

    wavefunc wf;

    // wavefunction coefficients
    psi::Psi<T>* psi = nullptr;

    std::string basisname; // PW or LCAO
    double esolver_KS_ne = 0.0;
    bool oscillate_esolver = false; // whether esolver is oscillated
#ifdef __MPI
    double iter_time;               // the start time of scf iteration
#else
    std::chrono::system_clock::time_point iter_time;
#endif
    double diag_ethr;               // the threshold for diagonalization
    double scf_thr;                 // scf density threshold
    double scf_ene_thr;             // scf energy threshold
    double drho;                    // the difference between rho_in (before HSolver) and rho_out (After HSolver)
    double hsolver_error;           // the error of HSolver
    int maxniter;                   // maximum iter steps for scf
    int niter;                      // iter steps actually used in scf
    int out_freq_elec;              // frequency for output
};
} // namespace ModuleESolver
#endif
