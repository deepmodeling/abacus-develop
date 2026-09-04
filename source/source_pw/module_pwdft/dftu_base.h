#ifndef DFTU_BASE_H
#define DFTU_BASE_H

#include "source_base/matrix.h"
#include "source_estate/occ_matrix.h"
#include "source_pw/module_pwdft/yukawa_screening.h"

#include <complex>
#include <memory>
#include <string>
#include <vector>


class UnitCell;
class Charge_Mixing;

class DFTUTest;

class Plus_U_Base
{
  friend class DFTUTest;

  public:
    Plus_U_Base();
    ~Plus_U_Base();

    /// allocate relevant data structures (base part, no LCAO types)
    void init_base(UnitCell& cell,
                   const int npol,
                   const int nspin,
                   const std::vector<int>& orbital_corr,
                   const bool yukawa_potential,
                   const double yukawa_lambda,
                   const std::string& global_readin_dir,
                   const std::string& global_out_dir,
                   const std::string& init_chg,
                   const std::string& device,
                   const int kpar,
                   const std::vector<double>& hubbard_u,
                   const double uramping,
                   const int occ_mat_ctrl,
                   const int mixing_dftu);

    void uramping_update();
    bool u_converged();

    // --- Accessors for U values and orbital configuration ---
    double get_u_current(int it) const { return u_current[it]; }
    int get_num_u_types() const { return static_cast<int>(u_current.size()); }
    int get_orbital_corr(int it) const { return orbital_corr[it]; }
    bool has_correlated_orbital(int it) const { return orbital_corr[it] != -1; }

    /// read-only access to the orbital_corr vector (length ntype)
    const std::vector<int>& get_orbital_corr_vec() const { return orbital_corr; }

    // --- Accessors for DFT+U configuration ---
    double get_uramping() const { return uramping; }
    int get_occ_mat_ctrl() const { return occ_mat_ctrl; }
    int get_cal_type() const { return cal_type; }
    bool use_yukawa() const { return yukawa_ != nullptr; }

    /// access the Yukawa screening object (non-null only when use_yukawa())
    YukawaScreening& yukawa() { return *yukawa_; }
    const YukawaScreening& yukawa() const { return *yukawa_; }

    void set_u_current(int it, double val) { u_current[it] = val; }

    double get_energy() const { return energy_u; }
    void set_energy(const double &e) { energy_u = e; }
    void set_double_energy() { energy_u *= 2.0; }

    /// interface for PW basis
    /// calculate the local occupation number matrix for PW based wave functions
    void cal_occ_pw(const void* psi_in,
                    const ModuleBase::matrix& wg_in,
                    const UnitCell& cell,
                    Charge_Mixing* p_chgmix,
                    const int* isk);

    /// get effective potential pointer for the given spin channel (PW basis)
    ///
    /// nspin=1: isk is ignored, returns &pot_uterm_pw[0]
    /// nspin=2: isk selects spin-up (0) or spin-down (1) half of the
    ///          split layout [all_up | all_dn]
    /// nspin=4: isk is ignored, returns &pot_uterm_pw[0] (all Pauli blocks)
    const std::complex<double>* get_pot_uterm_pw_spin(const int isk) const
    {
        if (nspin == 2 && isk == 1)
        {
            return pot_uterm_pw.data() + pot_uterm_pw.size() / 2;
        }
        return pot_uterm_pw.data();
    }

    /// get size of effective potential for a single spin channel (PW basis)
    ///
    /// nspin=1: full array size
    /// nspin=2: half of the total (one spin channel in split layout)
    /// nspin=4: full array size (all Pauli blocks are packed together)
    int get_size_pot_uterm_pw_spin() const
    {
        return (nspin == 2) ? static_cast<int>(pot_uterm_pw.size() / 2)
                            : static_cast<int>(pot_uterm_pw.size());
    }

    int get_size_pot_uterm_pw() const
    {
        return pot_uterm_pw.size();
    }

    // dftu can be calculated only after occ_mat has been initialized
    bool is_occ_mat_initialized() const { return occ_mat_initialized; }
    void mark_occ_mat_initialized() { occ_mat_initialized = true; }
    void mark_occ_mat_dirty() { occ_mat_initialized = false; }

    void enable_mixing() { mixing_dftu = 1; }

    /// direct access to the occupation matrix object (new write path)
    OccupationMatrix& occmat() { return occmat_; }
    const OccupationMatrix& occmat() const { return occmat_; }

  protected:
    // --- U values and orbital configuration (set in init_base) ---
    std::vector<double> u_current;
    std::vector<double> u_target;
    std::vector<int> orbital_corr;

    // --- DFT+U configuration flags ---
    double uramping = 0.0;
    int occ_mat_ctrl = 0;
    int mixing_dftu = 0;
    int nspin = 0;

    // --- State flags ---
    // dftu can be calculated only after occ_mat has been initialized
    bool occ_mat_initialized = false;

    // --- Occupation matrices ---
    OccupationMatrix occmat_;

    // --- Internal state ---
    double energy_u = 0.0;

    int cal_type = 3;
    std::string device;
    int kpar = 1;

    std::vector<std::complex<double>> pot_uterm_pw;
    std::vector<int> pot_uterm_pw_index;
    std::vector<double> uom_array;
    std::vector<double> uom_save;

    // Yukawa screening object; constructed only when use_yukawa() is true.
    // Owns the screening length, Slater integrals and derived U/J.
    std::unique_ptr<YukawaScreening> yukawa_;
};


#endif
