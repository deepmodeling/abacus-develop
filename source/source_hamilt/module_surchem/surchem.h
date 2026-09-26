#ifndef SURCHEM_H
#define SURCHEM_H

#include "source_base/atom_in.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "sccs/sccs_driver.h"

#include <iosfwd>

// forward-declared: used below only as pointer/reference
class Parallel_Grid;
class Structure_Factor;

/**
 * @brief Implicit-solvent settings, injected at the ESolver boundary.
 *
 * These deliberately carry no physical defaults. The only production instance of
 * surchem is ESolver_FP::solvent, which is always configured from Input_para via
 * surchem::set_parameters(); mirroring the INPUT defaults here would create a second
 * copy that could silently drift out of sync with input_parameter.h. Callers that
 * need specific values (unit tests included) must state them explicitly.
 */
struct SurchemParameters
{
    double eb_k = 0.0;    ///< relative permittivity of the bulk solvent
    double tau = 0.0;     ///< effective surface tension parameter
    double sigma_k = 0.0; ///< width of the diffuse cavity
    double nc_k = 0.0;    ///< cut-off charge density
    bool use_sccs = false;
    bool use_legacy_solvent = false;
    ModuleSccs::Boundary pcc_boundary = ModuleSccs::Boundary::Periodic;
    ModuleSccs::SccsConfig sccs_config;
    double expected_electron_count = 0.0;
    double expected_ionic_charge = 0.0;
    double normalization_tolerance = 1.0e-6;
    int pool_process_count = 1;
    double start_drho = 0.0;
    int start_nmax = 30;
    int debug = 0;
};

class surchem
{
  public:
    surchem();
    ~surchem();

    double* TOTN_real = nullptr;
    double* delta_phi = nullptr;
    double* epspot = nullptr;
    ModuleBase::matrix Vcav;
    ModuleBase::matrix Vel;
    double qs;

    static double Acav;
    static double Ael;

    // get atom info
    atom_in GetAtom;

    // allocate memory and deallocate them
    void allocate(const int& nrxx, const int& nspin);

    void clear();

    void set_parameters(const SurchemParameters& parameters);

    bool uses_sccs() const;
    bool uses_pcc() const;

    bool sccs_is_active() const;

    bool try_activate_sccs(int electronic_iteration, double drho);

    const ModuleSccs::SccsResult& sccs_result() const;

    void write_sccs_iteration(std::ostream& output) const;

    void write_sccs_diagnostics(std::ostream& output) const;

    void cal_epsilon(const ModulePW::PW_Basis* rho_basis, const double* PS_TOTN_real, double* epsilon, double* epsilon0);

    void cal_pseudo(const UnitCell& cell,
                    const Parallel_Grid& pgrid,
                    const ModulePW::PW_Basis* rho_basis,
                    const std::complex<double>* Porter_g,
                    std::complex<double>* PS_TOTN,
                    Structure_Factor* sf);

    void gauss_charge(const UnitCell& cell,
                      const Parallel_Grid& pgrid,
                      const ModulePW::PW_Basis* rho_basis,
                      std::complex<double>* N,
                      Structure_Factor* sf);

    void cal_totn(const UnitCell& cell,
                  const ModulePW::PW_Basis* rho_basis,
                  const std::complex<double>* Porter_g,
                  std::complex<double>* N,
                  std::complex<double>* TOTN,
                  const double* vlocal);

    void createcavity(const UnitCell& ucell,
                      const ModulePW::PW_Basis* rho_basis,
                      const std::complex<double>* PS_TOTN,
                      double* vwork);

    void cal_vcav(const UnitCell& ucell,
                  const ModulePW::PW_Basis* rho_basis,
                  std::complex<double>* PS_TOTN,
                  int nspin,
                  ModuleBase::matrix& v);

    void cal_vel(const UnitCell& cell,
                 const ModulePW::PW_Basis* rho_basis,
                 std::complex<double>* TOTN,
                 std::complex<double>* PS_TOTN,
                 int nspin,
                 ModuleBase::matrix& v);

    double cal_Ael(const UnitCell& cell,
                   const int& nrxx,  // num. of real space grids on current core
                   const int& nxyz); // total num. of real space grids

    double cal_Acav(const UnitCell& cell,
                    const int& nxyz); // total num. of real space grids

    void cal_Acomp(const UnitCell& cell,
                   const ModulePW::PW_Basis* rho_basis,
                   const double* const* const rho,
                   std::vector<double>& res);

    void minimize_cg(const UnitCell& ucell,
                     const ModulePW::PW_Basis* rho_basis,
                     double* d_eps,
                     const std::complex<double>* tot_N,
                     std::complex<double>* phi,
                     int& ncgsol);

    void Leps2(const UnitCell& ucell,
               const ModulePW::PW_Basis* rho_basis,
               std::complex<double>* phi,
               double* epsilon,            // epsilon from shapefunc, dim=nrxx
               std::complex<double>* gradphi_G_work,
               std::complex<double>* lp,
               ModuleBase::Vector3<double>* grad_phi_R,   // size: nrxx
               double* aux_R);

    void v_correction(const UnitCell& cell,
                      const Parallel_Grid& pgrid,
                      const ModulePW::PW_Basis* rho_basis,
                      const int& nspin,
                      const double* const* const rho,
                      const double* vlocal,
                      Structure_Factor* sf,
                      ModuleBase::matrix& v);

    void v_correction_sccs(const UnitCell& cell,
                           const ModulePW::PW_Basis& rho_basis,
                           int nspin,
                           const double* const* rho,
                           const double* vlocal,
                           ModuleBase::matrix& v);
    void v_correction_pcc(const UnitCell& cell,
                          const ModulePW::PW_Basis& rho_basis,
                          int nspin,
                          const double* const* rho,
                          ModuleBase::matrix& v);
    void cal_force_pcc(const UnitCell& cell, ModuleBase::matrix& force) const;

    void test_V_to_N(ModuleBase::matrix& v,
                     const UnitCell& cell,
                     const ModulePW::PW_Basis* rho_basis,
                     const double* const* const rho);

    void cal_force_sol(const UnitCell& cell,
                       const ModulePW::PW_Basis* rho_basis,
                       const ModuleBase::matrix& vloc,
                       int nspin,
                       ModuleBase::matrix& forcesol);

    void force_cor_one(const UnitCell& cell,
                       const ModulePW::PW_Basis* rho_basis,
                       const ModuleBase::matrix& vloc,
                       ModuleBase::matrix& forcesol);

    void force_cor_two(const UnitCell& cell,
                       const ModulePW::PW_Basis* rho_basis,
                       int nspin,
                       ModuleBase::matrix& forcesol);

    void cal_force_sccs(const UnitCell& cell,
                        const ModulePW::PW_Basis& rho_basis,
                        const ModuleBase::matrix& vloc,
                        ModuleBase::matrix& forcesol) const;

    void get_totn_reci(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, std::complex<double>* totn_reci);

    void induced_charge(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, double* induced_rho) const;

  private:
    SurchemParameters parameters_;
    bool parameters_set_ = false;
    bool sccs_active_ = false;
    ModuleSccs::SccsState sccs_state_;
    ModuleSccs::SccsResult sccs_result_;
    ModuleSccs::PccGeometry pcc_geometry_;
    ModuleSccs::Pcc2dGeometry pcc_2d_geometry_;
    ModuleSccs::MultipoleMoments pcc_moments_;
    ModuleSccs::Pcc2dMoments pcc_2d_moments_;
    ModuleSccs::Pcc2dMoments pcc_ionic_moments_2d_;
    double pcc_energy_rydberg_ = 0.0;
    double pcc_elapsed_seconds_ = 0.0;
    bool pcc_result_valid_ = false;
    double sccs_elapsed_seconds_ = 0.0;
};

#endif
