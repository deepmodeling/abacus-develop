#ifndef HSOLVER_H
#define HSOLVER_H

#include "diagh.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt/hamilt.h"
#include "module_psi/psi.h"
#include "src_pw/sto_wf.h"

#include <complex>

namespace hsolver
{

class HSolver
{
  public:
    HSolver(){};
    virtual ~HSolver(){
        delete pdiagh;
    };
    /*//initialization, used in construct function or restruct a new HSolver
    virtual void init(
        const Basis* pbas //We need Basis class here, use global class for this initialization first
        //const Input &in, //We need new Input class here, use global variable for this initialization first
        //elecstate::ElecState *pes
        )=0;
    //initialization, only be called for change some parameters only
    virtual void update(
        Input &in )=0;*/

    // solve Hamiltonian to electronic density in ElecState
    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<std::complex<double>>& ppsi, 
        elecstate::ElecState* pes, 
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }
    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<double>& ppsi, 
        elecstate::ElecState* pes, 
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }

    virtual void solve
    (
        hamilt::Hamilt* phm, 
        psi::Psi<std::complex<double>>& ppsi, 
        elecstate::ElecState* pes, 
        Stochastic_WF& stowf,
        const int istep,
        const int iter,
        const std::string method, 
        const bool skip_charge=false
    )
    {
        return;
    }

    std::string classname = "none";
    // choose method of DiagH for solve Hamiltonian matrix
    // cg, dav, elpa, scalapack, hpseps, cusolver
    std::string method = "none";
  public:
    double diag_ethr=0.0; //threshold for diagonalization
    //set diag_ethr according to drho
    //for lcao, we suppose the error is zero and we set diag_ethr to 0
    virtual double set_diagethr(const int istep, const int iter, const double drho)
    {
        return 0.0;
    }
    //reset diag_ethr according to drho and hsolver_error
    virtual double reset_diagethr(std::ofstream& ofs_running, const double hsover_error, const double drho)
    {
        return 0.0;
    }

    // calculate hsolver_error
    // for sdft and lcao, we suppose the error is zero 
    virtual double cal_hsolerror()
    {
        return 0.0;
    };

  protected:
    DiagH* pdiagh = nullptr; // for single Hamiltonian matrix diagonal solver

};

} // namespace hsolver
#endif