#ifndef DFTU_BASE_H
#define DFTU_BASE_H

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge_mixing.h"

#include <string>
#include <vector>


class DFTUTest;

class Plus_U_Base
{
  friend class DFTUTest;

  //=============================================================
  // public section
  //=============================================================
  public:
    Plus_U_Base();
    ~Plus_U_Base();

    /// allocate relevant data structures (base part, no LCAO types)
    void init_base(UnitCell& cell,
                   const int npol,
                   const int nspin,
                   const std::vector<int>& orbital_corr,
                   const bool yukawa_potential,
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
    double get_u_target(int it) const { return u_target[it]; }
    int get_num_u_types() const { return static_cast<int>(u_current.size()); }
    int get_orbital_corr(int it) const { return orbital_corr[it]; }
    bool has_correlated_orbital(int it) const { return orbital_corr[it] != -1; }
    const int* get_orbital_corr_data() const { return orbital_corr.data(); }

    /// read-only access to the orbital_corr vector (length ntype)
    const std::vector<int>& get_orbital_corr_vec() const { return orbital_corr; }

    /// read-only access to the iat->(l,n,m,ipol)->iwt lookup table
    const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>&
    get_iatlnmipol2iwt() const { return iatlnmipol2iwt; }

    // --- Accessors for DFT+U configuration ---
    double get_uramping() const { return uramping; }
    int get_occ_mat_ctrl() const { return occ_mat_ctrl; }
    int get_cal_type() const { return cal_type; }
    bool use_yukawa() const { return use_yukawa_; }

    double get_U_Yukawa(int it, int l, int n) const { return U_Yukawa[it][l][n]; }
    double get_J_Yukawa(int it, int l, int n) const { return J_Yukawa[it][l][n]; }
    void set_U_Yukawa(int it, int l, int n, double val) { U_Yukawa[it][l][n] = val; }
    void set_J_Yukawa(int it, int l, int n, double val) { J_Yukawa[it][l][n] = val; }
    double get_lambda() const { return lambda; }
    void set_lambda(double l) { lambda = l; }
    std::vector<std::vector<std::vector<std::vector<double>>>>& get_Fk_data() { return Fk; }
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

    /// get effective potential matrix for PW base (per-atom, raw index)
    /// @deprecated Use get_pot_uterm_pw_spin() for nspin-aware access.
    [[deprecated("Use get_pot_uterm_pw_spin() for nspin-aware access")]]
    const std::complex<double>* get_pot_uterm_pw(const int iat) const
    {
        return &(pot_uterm_pw[pot_uterm_pw_index[iat]]);
    }

    int get_size_pot_uterm_pw() const
    {
        return pot_uterm_pw.size();
    }

    // dftu can be calculated only after occ_mat has been initialized
    bool is_occ_mat_initialized() const { return occ_mat_initialized; }
    void mark_occ_mat_initialized() { occ_mat_initialized = true; }
    void mark_occ_mat_dirty() { occ_mat_initialized = false; }

    bool is_mixing_enabled() const { return mixing_dftu != 0; }
    void enable_mixing() { mixing_dftu = 1; }

    /// get occupation matrix element occ_mat[iat][l][n][spin](m1,m2)
    double get_occ_mat(const int iat, const int l, const int n, const int spin,
                     const int m1, const int m2) const
    {
        return occ_mat[iat][l][n][spin](m1, m2);
    }

    /// get saved occupation matrix element occ_mat_save[iat][l][n][spin](m1,m2)
    double get_occ_mat_save(const int iat, const int l, const int n, const int spin,
                          const int m1, const int m2) const
    {
        return occ_mat_save[iat][l][n][spin](m1, m2);
    }

    /// set occupation matrix element occ_mat[iat][l][n][spin](m1,m2)
    void set_occ_mat(const int iat, const int l, const int n, const int spin,
                     const int m1, const int m2, const double val)
    {
        occ_mat[iat][l][n][spin](m1, m2) = val;
    }

    /// get reference to occ_mat data
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& get_occ_mat_data() { return occ_mat; }
    /// get reference to occ_mat_save data
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& get_occ_mat_save_data() { return occ_mat_save; }
    /// get occ_mat_initialized flag
    bool get_occ_mat_initialized() const { return occ_mat_initialized; }
    /// set occ_mat_initialized flag
    void set_occ_mat_initialized(bool val) { occ_mat_initialized = val; }

    /// get flat occupation matrix for an atom's correlated orbital.
    /// nspin=1: fills occ with occ_mat[iat][l][0][0] data
    /// nspin=2: fills occ with interleaved occ_mat[iat][l][0][0] and [1] data
    /// nspin=4: fills occ with occ_mat[iat][l][0][0] data (all 4 Pauli blocks)
    void get_occ_mat_flat(const int iat, const int l, std::vector<double>& occ) const;

    /// set flat occupation matrix for an atom's correlated orbital (write-back)
    void set_occ_mat_flat(const int iat, const int l, const int spin,
                        const std::vector<double>& occ);

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
    bool use_yukawa_ = false;

    // --- State flags ---
    // dftu can be calculated only after occ_mat has been initialized
    bool occ_mat_initialized = false;

    // --- Occupation matrices ---
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> occ_mat;
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> occ_mat_save;

    // --- Internal state ---
    double energy_u = 0.0;

    int cal_type = 3;
    std::string device;
    int kpar = 1;

    // transform between iwt index and it, ia, L, N and m index
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>
        iatlnmipol2iwt;

    void copy_occ_mat(const UnitCell& ucell);
    void zero_occ_mat(const UnitCell& ucell);
    void mix_occ_mat(const UnitCell& ucell, const double& mixing_beta);
    void set_occ_mat(const UnitCell& ucell);

    /// accumulate occ_mat from psi for all k-points (per-device template)
    template <typename Device>
    void accumulate_occ_one_k(const void* psi_in,
                              const ModuleBase::matrix& wg_in,
                              const UnitCell& cell,
                              const int* isk);

    /// reduce occ_mat across k-pools (per-atom, nspin-aware)
    void reduce_occ_mat(const UnitCell& cell);

    /// copy occ_mat to uom_array for mixing (nspin-aware split layout)
    void sync_occ_to_uom(const UnitCell& cell);

    /// compute effective potential pot_onsite and DFT+U energy from occ_mat
    /// (assumes occ_mat has already been reduced across k-pools)
    void compute_eff_pot_and_energy(const UnitCell& cell);

    std::vector<std::complex<double>> pot_uterm_pw;
    std::vector<int> pot_uterm_pw_index;
    std::vector<double> uom_array;
    std::vector<double> uom_save;

    // Yukawa-related members (base part, no LCAO dependency)
    double lambda = 0.0;
    std::vector<std::vector<std::vector<std::vector<double>>>> Fk;
    std::vector<std::vector<std::vector<double>>> U_Yukawa;
    std::vector<std::vector<std::vector<double>>> J_Yukawa;

    void read_occup_m(const UnitCell& ucell,
                      const std::string& fn,
                      const std::string& init_chg,
                      int nspin,
                      int npol);
    void local_occup_bcast(const UnitCell& ucell,
                           int nspin,
                           int npol);
};


#endif
