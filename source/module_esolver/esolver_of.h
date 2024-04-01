#ifndef ESOLVER_OF_H
#define ESOLVER_OF_H

#include "esolver_fp.h"
#include "module_base/opt_DCsrch.h"
#include "module_base/opt_TN.hpp"
#include "module_elecstate/module_charge/charge.h"
#include "module_elecstate/module_charge/charge_extra.h" // liuyu add 2022-11-07
#include "module_hamilt_pw/hamilt_ofdft/kedf_lkt.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_tf.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_vw.h"
#include "module_hamilt_pw/hamilt_ofdft/kedf_wt.h"
#include "module_psi/psi.h"

namespace ModuleESolver
{
class ESolver_OF : public ESolver_FP
{
  public:
    ESolver_OF();
    ~ESolver_OF();

    virtual void init(Input& inp, UnitCell& ucell) override;

    virtual void init_after_vc(Input& inp, UnitCell& ucell) override;

    virtual void run(int istep, UnitCell& ucell) override;

    virtual void post_process() override;

    virtual double cal_energy() override;

    virtual void cal_force(ModuleBase::matrix& force) override;

    virtual void cal_stress(ModuleBase::matrix& stress) override;

    virtual int getniter() override
    {
        return this->iter_;
    }

  private:
    // ======================= variables ==========================
    // ---------- the kinetic energy density functionals ----------
    KEDF_TF* tf_ = nullptr;
    KEDF_vW* vw_ = nullptr;
    KEDF_WT* wt_ = nullptr;
    KEDF_LKT* lkt_ = nullptr;

    // charge extrapolation liuyu add 2022-11-07
    Charge_Extra CE_;

    // ----------------- the optimization methods ------------------
    ModuleBase::Opt_CG* opt_cg_ = nullptr;
    ModuleBase::Opt_TN* opt_tn_ = nullptr;
    ModuleBase::Opt_DCsrch* opt_dcsrch_ = nullptr;
    ModuleBase::Opt_CG* opt_cg_mag_ = nullptr; // for spin2 case, under testing

    // ----------------- necessary parameters from INPUT ------------
    std::string of_kinetic_ = "wt";  // Kinetic energy functional, such as TF, VW, WT
    std::string of_method_ = "tn";   // optimization method, include cg1, cg2, tn (default), bfgs
    std::string of_conv_ = "energy"; // select the convergence criterion, potential, energy (default), or both
    double of_tole_ = 2e-6; // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp_ = 1e-5; // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    int max_iter_ = 50;     // scf_nmax

    // ------------------ parameters from other module --------------
    double dV_ = 0;           // volume of one grid point in real space
    double* nelec_ = nullptr; // number of electrons with each spin

    // ----- parameters and arrays used in density optimization -----
    int iter_ = 0;                                // iteration number
    double** pdirect_ = nullptr;                  // optimization direction of phi, which is sqrt(rho)
    std::complex<double>** precip_dir_ = nullptr; // direction in reciprocal space, used when of_full_pw=false.
    double* theta_ = nullptr;                     // step length
    double** pdEdphi_ = nullptr;                  // dE/dphi
    double** pdLdphi_ = nullptr;                  // dL/dphi
    double** pphi_ = nullptr;                     // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    char* task_ = nullptr;                        // used in line search
    double* mu_ = nullptr;                        // chemical potential
    int tn_spin_flag_ = -1;                       // spin flag used in cal_potential, which will be called by opt_tn
    int max_dcsrch_ = 200;                        // max no. of line search
    int flag_ = -1;                               // flag of TN
    Charge* ptemp_rho_ = nullptr;                 // used in line search
    psi::Psi<double>* psi_ = nullptr;             // sqrt(rho)

    // ----------------- used for convergence check -------------------
    bool conv_ = false;
    double energy_llast_ = 0;
    double energy_last_ = 0;
    double energy_current_ = 0;
    double normdLdphi_llast_ = 100;
    double normdLdphi_last_ = 100;
    double normdLdphi_ = 100.;

    // ==================== main process of OFDFT ======================
    void before_opt(const int istep, UnitCell& ucell);
    void update_potential(UnitCell& ucell);
    void optimize(UnitCell& ucell);
    void update_rho();
    bool check_exit();
    void after_opt(const int istep, UnitCell& ucell);

    // ============================ tools ===============================
    // --------------------- initialize ---------------------------------
    void init_elecstate(UnitCell& ucell);
    void allocate_array();

    // --------------------- calculate physical qualities ---------------
    void cal_potential(double* ptemp_phi, double* rdLdphi);
    void cal_dEdtheta(double** ptemp_phi, Charge* temp_rho, UnitCell& ucell, double* ptheta, double* rdEdtheta);
    double cal_mu(double* pphi, double* pdEdphi, double nelec);

    // --------------------- determine the optimization direction -------
    void adjust_direction();
    void check_direction(double* dEdtheta, double** ptemp_phi, UnitCell& ucell);
    void test_direction(double* dEdtheta, double** ptemp_phi, UnitCell& ucell);

    // --------------------- output the necessary information -----------
    void print_info();

    // --------------------- interface to blas --------------------------
    double inner_product(double* pa, double* pb, int length, double dV = 1)
    {
        double innerproduct = BlasConnector::dot(length, pa, 1, pb, 1);
        innerproduct *= dV;
        return innerproduct;
    }

    // ---------------------- interfaces to KEDF ------------------------
    void init_kedf(Input& inp);
    void kinetic_potential(double** prho, double** pphi, ModuleBase::matrix& rpot);
    double kinetic_energy();
    void kinetic_stress(ModuleBase::matrix& kinetic_stress);

    // ---------------------- interfaces to optimization methods --------
    void init_opt();
    void get_direction();
    void get_step_length(double* dEdtheta, double** ptemp_phi, UnitCell& ucell);
};
} // namespace ModuleESolver

#endif
