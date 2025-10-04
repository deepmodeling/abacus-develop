#ifndef ESOLVER_KS_H
#define ESOLVER_KS_H

#include <cstring>
//#include <fstream>

// for first-principles esolver
#include "esolver_fp.h"
// for plane wave basis set 
#include "source_basis/module_pw/pw_basis_k.h"
// for k-points in Brillouin zone
#include "source_cell/klist.h"
// for charge mixing
#include "source_estate/module_charge/charge_mixing.h"
// for electronic wave functions
#include "source_psi/psi.h"
// for Hamiltonian
#include "source_hamilt/hamilt.h"

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

    virtual void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    virtual void runner(UnitCell& ucell, const int istep) override;

    virtual void after_all_runners(UnitCell& ucell) override;

  protected:
    //! Something to do before SCF iterations.
    virtual void before_scf(UnitCell& ucell, const int istep) override;

    //! Something to do before hamilt2rho function in each iter loop.
    virtual void iter_init(UnitCell& ucell, const int istep, const int iter);

    //! Something to do after hamilt2rho function in each iter loop.
    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    // calculate electron density from a specific Hamiltonian with ethr
    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr);

    // calculate electron density from a specific Hamiltonian
    void hamilt2rho(UnitCell& ucell, const int istep, const int iter, const double ethr);

    //! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    //! <Temporary> It should be replaced by a function in Hamilt Class
    virtual void update_pot(UnitCell& ucell, const int istep, const int iter, const bool conv_esolver){};

    //! Hamiltonian
    hamilt::Hamilt<T, Device>* p_hamilt = nullptr;

    //! PW for wave functions, only used in KSDFT, not in OFDFT
    ModulePW::PW_Basis_K* pw_wfc = nullptr;

    //! Charge mixing method
    Charge_Mixing* p_chgmix = nullptr;

    //! nonlocal pseudopotentials
    pseudopot_cell_vnl ppcell;

    //! Electronic wavefunctions
    psi::Psi<T>* psi = nullptr;

    std::string basisname;      //! esolver_ks_lcao.cpp
    double esolver_KS_ne = 0.0; //! number of electrons
    double diag_ethr;           //! the threshold for diagonalization
    double scf_thr;             //! scf density threshold
    double scf_ene_thr;         //! scf energy threshold
    double drho;                //! the difference between rho_in (before HSolver) and rho_out (After HSolver)
    double hsolver_error;       //! the error of HSolver
    int maxniter;               //! maximum iter steps for scf
    int niter;                  //! iter steps actually used in scf
    bool oscillate_esolver = false; // whether esolver is oscillated

    bool scf_nmax_flag = false; // whether scf has reached nmax, mohan add 20250921
};
} // namespace ModuleESolver
#endif
