#ifndef SYMM_RHO_H
#define SYMM_RHO_H
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_estate/module_charge/charge.h"
#include "source_pw/module_pwdft/parallel_grid.h"

class Symmetry_rho
{
  public:
    Symmetry_rho();
    ~Symmetry_rho();

    /**
     * @brief Symmetrize charge density for all spin channels
     * 
     * This is a static helper function that symmetrizes the charge density
     * for all spin channels by calling begin() for each spin.
     * 
     * @param nspin Number of spin channels
     * @param chr Charge object containing the density
     * @param pw Plane wave basis
     * @param symm Symmetry object
     */
    static void symmetrize_rho(const int nspin,
                               const Charge& chr,
                               const ModulePW::PW_Basis* pw,
                               ModuleSymmetry::Symmetry& symm);

    void begin(const int& spin_now,
               const Charge& CHR,
               const ModulePW::PW_Basis* pw,
               ModuleSymmetry::Symmetry& symm) const;

    void begin(const int& spin_now,
               double** rho,
               std::complex<double>** rhog,
               int ngmc,
               double** kin_r,
               const ModulePW::PW_Basis* pw,
               ModuleSymmetry::Symmetry& symm) const;

    /// @brief Symmetrize the nspin=4 spin density (rho^x, rho^y, rho^z = rho[1,2,3]) with the
    ///        coupled spin rotation. The charge component rho^0 = rho[0] is handled separately
    ///        by the ordinary scalar begin().
    void begin_soc(const Charge& CHR,
                   const ModulePW::PW_Basis* pw,
                   ModuleSymmetry::Symmetry& symm) const;

  private:
    // in real space:
    void psymm(double* rho_part,
               const ModulePW::PW_Basis* pw,
               Parallel_Grid& Pgrid,
               ModuleSymmetry::Symmetry& symm) const;
    // in reciprocal space:
    void psymmg(std::complex<double>* rhog_part,
                const ModulePW::PW_Basis* rho_basis,
                ModuleSymmetry::Symmetry& symm) const;
    // in reciprocal space, the three coupled spin components (rho^x, rho^y, rho^z) for nspin=4:
    void psymmg_soc(std::complex<double>* rhog_x,
                    std::complex<double>* rhog_y,
                    std::complex<double>* rhog_z,
                    const ModulePW::PW_Basis* rho_basis,
                    ModuleSymmetry::Symmetry& symm) const;
#ifdef __MPI
    void reduce_to_fullrhog(const ModulePW::PW_Basis* rho_basis,
                            std::complex<double>* rhogtot,
                            std::complex<double>* rhogin,
                            int* ig2isztot,
                            const int* ig2iszin,
                            int max_npw) const;
    void rhog_piece_to_all(const ModulePW::PW_Basis* rho_basis,
                           std::complex<double>* rhogtot,
                           std::complex<double>* rhog_part) const;
#endif
    void get_ixyz2ipw(const ModulePW::PW_Basis* rho_basis,
                      const int* ig2isztot,
                      const int* fftixy2is,
                      int* ixyz2ipw) const; //(ix, iy, iz) -> (ip, ig)
};

#endif
