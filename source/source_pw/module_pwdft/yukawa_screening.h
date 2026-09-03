#ifndef YUKAWA_SCREENING_H
#define YUKAWA_SCREENING_H

#include <vector>

class UnitCell;
class LCAO_Orbitals;

/**
 * @brief Yukawa-screened DFT+U: self-consistent U/J from Slater integrals.
 *
 * Encapsulates the Yukawa screening length (lambda), the Slater integrals Fk,
 * and the derived U_Yukawa / J_Yukawa values. All Yukawa-related state that
 * used to live on Plus_U_Base / Plus_U is owned here, so the DFT+U classes
 * only hold an instance of this class when Yukawa screening is enabled.
 *
 * Currently only the LCAO path drives the calculation (it needs the radial
 * orbitals from LCAO_Orbitals), but the class itself has no LCAO-only data.
 */
class YukawaScreening
{
  public:
    YukawaScreening() = default;
    ~YukawaScreening() = default;

    /// allocate Fk / U_Yukawa / J_Yukawa according to the cell and record the
    /// user-provided screening length (yukawa_lambda_cfg > 0 means fixed).
    void init(const UnitCell& cell,
              const std::vector<int>& orbital_corr,
              double yukawa_lambda_cfg);

    /// determine lambda: use the fixed config value when positive, otherwise
    /// estimate from the charge density (Thomas-Fermi-like) and rescale by 1.6.
    void cal_lambda(double** rho, int nrxx, int nspin);

    /// compute Slater integrals Fk for the correlated orbital of atom type T.
    void cal_slater_Fk(const UnitCell& ucell, int L, int T, const LCAO_Orbitals* orb);

    /// drive cal_lambda + cal_slater_Fk over all correlated orbitals and derive
    /// U_Yukawa / J_Yukawa. Returns via get_U/get_J; u_current of the owning
    /// DFT+U object is updated by the caller.
    void cal_slater_UJ(const UnitCell& ucell,
                       double** rho,
                       int nrxx,
                       int nspin,
                       const LCAO_Orbitals* orb);

    double get_lambda() const { return lambda_; }
    double get_U(int it, int l, int n) const { return U_Yukawa_[it][l][n]; }
    double get_J(int it, int l, int n) const { return J_Yukawa_[it][l][n]; }
    /// effective U-J of the correlated orbital (n = 0) for atom type it
    double get_Ueff(int it) const
    {
        const int l = orbital_corr_[it];
        return U_Yukawa_[it][l][0] - J_Yukawa_[it][l][0];
    }

  private:
    /// spherical modified Bessel function of the first kind, orders 0/2/4/6
    static double spherical_Bessel(int k, double r, double lambda);
    /// spherical modified Hankel function of the second kind, orders 0/2/4/6
    static double spherical_Hankel(int k, double r, double lambda);

    double lambda_ = 0.0;
    double yukawa_lambda_cfg_ = 0.0;
    std::vector<int> orbital_corr_;
    std::vector<std::vector<std::vector<std::vector<double>>>> Fk_;
    std::vector<std::vector<std::vector<double>>> U_Yukawa_;
    std::vector<std::vector<std::vector<double>>> J_Yukawa_;
};

#endif
