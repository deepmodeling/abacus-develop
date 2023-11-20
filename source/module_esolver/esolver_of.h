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
    ESolver_OF()
    {
        this->classname = "ESolver_OF";
        this->task_ = new char[60];
    }

    ~ESolver_OF()
    {
        delete psi_;
        delete[] this->pphi_;

        for (int i = 0; i < GlobalV::NSPIN; ++i)
        {
            delete[] this->pdirect_[i];
            delete[] this->pdLdphi_[i];
            delete[] this->pdEdphi_[i];
            delete[] this->precipDir_[i];
        }
        delete[] this->pdirect_;
        delete[] this->pdLdphi_;
        delete[] this->pdEdphi_;
        delete[] this->precipDir_;

        delete[] this->nelec_;
        delete[] this->theta_;
        delete[] this->mu_;
        delete[] this->task_;
        delete this->opt_cg_mag;
        delete this->ptempRho_;

        delete this->tf_;
        delete this->vw_;
        delete this->wt_;
        delete this->lkt_;
    }

    virtual void Init(Input &inp, UnitCell &ucell) override;
    virtual void init_after_vc(Input &inp, UnitCell &ucell) override;
    virtual void Run(int istep, UnitCell& ucell) override;
    virtual void postprocess() override;

    virtual double cal_Energy() override;
    virtual void cal_Force(ModuleBase::matrix &force) override;
    virtual void cal_Stress(ModuleBase::matrix &stress) override;

    virtual int getniter() override {
        return this->iter_;
    }

private:
    // kinetic energy density functionals
    KEDF_TF* tf_ = nullptr;
    KEDF_vW* vw_ = nullptr;
    KEDF_WT* wt_ = nullptr;
    KEDF_LKT* lkt_ = nullptr;

    // charge extrapolation liuyu add 2022-11-07
    Charge_Extra CE;
    psi::Psi<double>* psi_=nullptr;

    // optimization methods
    ModuleBase::Opt_CG opt_cg;
    ModuleBase::Opt_TN opt_tn;
    ModuleBase::Opt_DCsrch opt_dcsrch;
    ModuleBase::Opt_CG *opt_cg_mag = nullptr; // for spin2 case, under testing

    // from Input
    std::string of_kinetic_ = "wt";   // Kinetic energy functional, such as TF, VW, WT
    std::string of_method_ = "tn";    // optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv_ = "energy";  // select the convergence criterion, potential, energy (default), or both
    double of_tole_ = 2e-6;      // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp_ = 1e-5;      // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    int max_iter_ = 50;           // scf_nmax

    // parameters from other module
    double dV_ = 0; // CELL
    double *nelec_ = nullptr;              // number of electrons with each spin

    // used in density optimization
    int iter_ = 0;                               // iteration number
    double **pdirect_ = nullptr;                    // optimization direction of phi, which is sqrt(rho)
    std::complex<double> **precipDir_ = nullptr;    // direction in reciprocal space, used when of_full_pw=false.
    double *theta_ = nullptr;                       // step length
    double **pdEdphi_ = nullptr;                    // dE/dphi
    double **pdLdphi_ = nullptr;                    // dL/dphi
    double **pphi_ = nullptr;                       // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    char *task_ = nullptr;                          // used in line search
    double *mu_ = nullptr;                          // chemical potential
    int tnSpinFlag_ = -1;                        // spin flag used in calV, which will be called by opt_tn
    int maxDCsrch_ = 200;                        // max no. of line search
    int flag_ = -1;                              // flag of TN

    Charge* ptempRho_ = nullptr;                 // used in line search

    // // test rho convergence criterion
    // double *pdeltaRhoHar = nullptr; // 4pi*rhog/k^2
    // double deltaRhoG = 0.; // 1/2\iint{deltaRho(r)deltaRho(r')/|r-r'|drdr'}
    // double deltaRhoR = 0.; // \int{|deltaRho(r)|dr}

    // used in convergence check
    bool conv_ = false;
    double energy_llast_ = 0;
    double energy_last_ = 0;
    double energy_current_ = 0;
    double normdLdphi_llast_ = 100;
    double normdLdphi_last_ = 100;
    double normdLdphi_ = 100.;

    // main process of OFDFT
    void beforeOpt(const int istep);
    void updateV();
    void solveV();
    void getNextDirect();
    void updateRho();
    bool checkExit();
    void printInfo();
    void afterOpt(const int istep);

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
    void init_kedf();
    void kineticPotential(double **prho, double **pphi, ModuleBase::matrix &rpot);
    double kineticEnergy();
    void kinetic_stress(ModuleBase::matrix &kinetic_stress);
};
}

#endif