#ifndef ESOLVER_OF_H
#define ESOLVER_OF_H

#include "esolver_fp.h"
#include "module_base/opt_TN.hpp"
#include "module_base/opt_DCsrch.h"
#include "module_psi/psi.h"
#include "module_elecstate/module_charge/charge_extra.h"    // liuyu add 2022-11-07
#include "module_hamilt_pw/hamilt_ofdft/kedf_tf.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_vw.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_wt.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_lkt.h"
#include "module_elecstate/module_charge/charge.h"

namespace ModuleESolver
{
class ESolver_OF: public ESolver_FP
{
// =========== TO DO LIST =============
// MORE KEDF
// GAMMA ONLY
// SPIN POLARISE
public:
    psi::Psi<double>* psi=nullptr;

    ESolver_OF()
    {
        this->classname = "ESolver_OF";
        this->task = new char[60];
    }

    ~ESolver_OF()
    {
        if(this->psi != nullptr)
        {
            delete psi;
        }
        if (this->pdirect != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdirect[i];
            }
            delete[] this->pdirect;
        } 
        if (this->pphi != NULL) delete[] this->pphi;
        if (this->pdLdphi != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdLdphi[i];
            }
            delete[] this->pdLdphi;
        } 
        if (this->pdEdphi != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdEdphi[i];
            }
            delete[] this->pdEdphi;
        } 
        if (this->precipDir != NULL)
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                delete[] this->precipDir[is];
            }
            delete[] this->precipDir;
        }

        if (this->nelec != NULL) delete[] this->nelec;
        if (this->theta != NULL) delete[] this->theta;
        if (this->mu != NULL) delete[] this->mu;
        if (this->task != NULL) delete[] this->task;
        if (this->opt_cg_mag != NULL) delete this->opt_cg_mag;
        delete this->ptempRho;
    }

    virtual void Init(Input &inp, UnitCell &ucell) override;
    virtual void init_after_vc(Input &inp, UnitCell &ucell) override;
    virtual void Run(int istep, UnitCell& ucell) override;
    virtual void postprocess() override;

    virtual double cal_Energy() override;
    virtual void cal_Force(ModuleBase::matrix &force) override;
    virtual void cal_Stress(ModuleBase::matrix &stress) override;

    virtual int getniter() override {
        return this->iter;
    }

private:
    // kinetic energy density functionals
    KEDF_TF tf;
    KEDF_vW vw;
    KEDF_WT wt;
    KEDF_LKT lkt;

    // charge extrapolation liuyu add 2022-11-07
    Charge_Extra CE;

    // optimization methods
    ModuleBase::Opt_CG opt_cg;
    ModuleBase::Opt_TN opt_tn;
    ModuleBase::Opt_DCsrch opt_dcsrch;
    ModuleBase::Opt_CG *opt_cg_mag = NULL; // for spin2 case, under testing

    // from Input
    std::string of_kinetic = "wt";   // Kinetic energy functional, such as TF, VW, WT
    std::string of_method = "tn";    // optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv = "energy";  // select the convergence criterion, potential, energy (default), or both
    double of_tole = 2e-6;      // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp = 1e-5;      // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    int maxIter = 50;           // scf_nmax

    // parameters from other module
    int nrxx = 0; // PWBASIS
    double dV = 0; // CELL
    double *nelec = NULL;              // number of electrons with each spin

    // used in density optimization
    int iter = 0;                               // iteration number
    double **pdirect = NULL;                    // optimization direction of phi, which is sqrt(rho)
    std::complex<double> **precipDir = NULL;    // direction in reciprocal space, used when of_full_pw=false.
    double *theta = NULL;                       // step length
    double **pdEdphi = NULL;                    // dE/dphi
    double **pdLdphi = NULL;                    // dL/dphi
    double **pphi = NULL;                       // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    char *task = NULL;                          // used in line search
    double *mu = NULL;                          // chemical potential
    int tnSpinFlag = -1;                        // spin flag used in calV, which will be called by opt_tn
    int maxDCsrch = 200;                        // max no. of line search
    int flag = -1;                              // flag of TN

    Charge* ptempRho = nullptr;                 // used in line search

    // // test rho convergence criterion
    // double *pdeltaRhoHar = NULL; // 4pi*rhog/k^2
    // double deltaRhoG = 0.; // 1/2\iint{deltaRho(r)deltaRho(r')/|r-r'|drdr'}
    // double deltaRhoR = 0.; // \int{|deltaRho(r)|dr}

    // used in convergence check
    bool conv = false;
    double energy_llast = 0;
    double energy_last = 0;
    double energy_current = 0;
    double normdLdphi_llast = 100;
    double normdLdphi_last = 100;
    double normdLdphi = 100.;

    // main process of OFDFT
    void beforeOpt(const int istep);
    void updateV();
    void solveV();
    void getNextDirect();
    void updateRho();
    bool checkExit();
    void printInfo();
    void afterOpt();

    // tools
    void calV(double *ptempPhi, double *rdLdphi);
    void caldEdtheta(double **ptempPhi, Charge* ptempRho, double *ptheta, double *rdEdtheta);
    double cal_mu(double *pphi, double *pdEdphi, double nelec);
    double inner_product(double *pa, double *pb, int length, double dV=1)
    {
        double innerproduct = 0.;
        for (int i = 0; i < length; ++i) innerproduct += pa[i] * pb[i];
        innerproduct *= dV;
        return innerproduct;
    }

    // interfaces to KEDF
    void kineticPotential(double **prho, double **pphi, ModuleBase::matrix &rpot);
    double kineticEnergy();
};
}

#endif