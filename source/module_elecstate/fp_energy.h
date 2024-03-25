/**
 * @file fp_energy.h
 * @brief This file contains all energies about first-principle calculations
 */
#ifndef FP_ENERGY_H
#define FP_ENERGY_H
namespace elecstate
{
/**
 * @struct fp_energy
 * @brief energies contribute to Hamholtz free energy
 */
struct fenergy
{
    double etot = 0.0;     ///< the total free energy
    double etot_old = 0.0; ///< old total free energy
    double etot_delta = 0.0; // the difference of total energy between two steps = etot - etot_old

    double eband = 0.0;          ///< the band energy
    double deband = 0.0;         ///< correction for variational energy
    double etxc = 0.0;           ///< E_xc[n(r)] exchange and correlation energy
    double etxcc = 0.0;          ///< the nlcc exchange and correlation
    double vtxc = 0.0;           ///< v_xc(r) = \frac{\delta E_xc}{\delta n} exchange and correlation potential
    double ewald_energy = 0.0;   ///< Ewald energy (ion-ion)
    double hartree_energy = 0.0; ///< Hartree energy (electron-electron)
    double demet = 0.0;          ///< correction for metals or entropy (-TS)
    double descf = 0.0;          ///< correction by the difference of rho
    double exx = 0.0;            ///< the exact exchange energy.

    double efield = 0.0;    ///< dipole potential in surface calculations
    double gatefield = 0.0; ///< correction energy for gatefield
    double evdw = 0.0;      ///< the vdw energy

    double etot_harris = 0.0;   ///< total energy of harris functional
    double deband_harris = 0.0; ///< correction for harris energy

    double esol_el = 0.0;  ///< the implicit solvation energy Ael
    double esol_cav = 0.0; ///< the implicit solvation energy Acav

    double edftu = 0.0;       ///< DFT+U energy
    double edeepks_scf = 0.0; /// DeePKS energy

    double escon = 0.0; ///< spin constraint energy

    double ekinetic = 0.0;      /// kinetic energy, used in OFDFT
    double eion_elec = 0.0;     /// ion-electron interaction energy, used in OFDFT

    double calculate_etot();
    double calculate_harris();
    void clear_all();
    void print_all() const;
};

/**
 * @struct efermi
 * @brief Fermi energies
 */
struct efermi
{
    double ef = 0.0;         ///< Fermi energy
    double ef_up = 0.0;      ///< spin up Fermi energy
    double ef_dw = 0.0;      ///< spin down Fermi energy
    bool two_efermi = false; ///<
    double& get_ef(const int& is);
    double get_efval(const int& is) const;
};

} // namespace elecstate
#endif