#ifndef ESOLVER_SDFT_PW_H
#define ESOLVER_SDFT_PW_H

#include "esolver_ks_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
#include "module_hamilt_pw/hamilt_stodft/sto_hchi.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"

namespace ModuleESolver
{

class ESolver_SDFT_PW: public ESolver_KS_PW<double>
{
public:
    ESolver_SDFT_PW();
    ~ESolver_SDFT_PW();
    void Init(Input &inp, UnitCell &cell) override;
    double cal_Energy() override;
    void cal_Force(ModuleBase::matrix& force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;
public:
    Stochastic_WF stowf;
    Stochastic_Iter stoiter;

protected:
    virtual void beforescf(const int istep) override; 
    // virtual void eachiterinit(int iter) override; 
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void nscf() override;
    virtual void othercalculation(const int istep) override;
    virtual void eachiterfinish(const int iter) override; 
    virtual void afterscf(const int istep) override;
    virtual void postprocess() override;

public:
    /**
     * @brief calculate Stochastic Kubo-Greenwood
     * 
     * @param nche_KG Number Chebyshev orders
     * @param fwhmin FWHM
     * @param wcut cutoff omega
     * @param dw_in omega step
     * @param dt_in t step
     * @param nbatch t step batch
     * @param npart_sto number stochastic wavefunctions parts to evalution at the same time  
     */
    void sKG(const int nche_KG, const double fwhmin, const double wcut, 
             const double dw_in, const double dt_in, const int nbatch, const int npart_sto);
    //calculate DOS
    void caldos(const int nche_dos, const double sigmain, 
            const double emin, const double emax, const double de, const int npart);

private:
    int nche_sto; ///<norder of Chebyshev
    int method_sto; ///<method of SDFT
    
    /**
     * @brief Check if Emin and Emax are converged
     * 
     * @param nche_in N order of Chebyshev expansion
     */
    void check_che(const int nche_in);

    /**
     * @brief Set the N order of Chebyshev expansion for conductivities
     * 
     * @param dt t step
     * @param nbatch number of t batch
     * @param cond_thr threshold of errors for conductivities
     * @return N order of Chebyshev
     */
    int set_cond_nche(const double dt, const int nbatch, const double cond_thr);

    /**
     * @brief calculate thermal current and electrical current
     * 
     */
    void cal_j(const psi::Psi<std::complex<double>>& psi_in,
               psi::Psi<std::complex<double>>& j1psi,
               psi::Psi<std::complex<double>>& j2psi,
               hamilt::Velocity& velop,
               const double& mu,
               const int& start_band,
               const int& nbands,
               const int& npw);


};

}//namespace ModuleESolver

//temporary setting: removed GlobalC but not breaking design philosophy
namespace GlobalTemp
{

extern const ModuleBase::matrix* veff;

}

#endif