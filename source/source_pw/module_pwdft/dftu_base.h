#ifndef DFTU_BASE_H
#define DFTU_BASE_H

#include "source_base/matrix.h"
#include "source_estate/occ_matrix.h"
#include "source_estate/occ_mixer.h"
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
    // DFT+U formalism & double-counting (DC) scheme.
    // Only dud_fll is implemented; the others are placeholders (return 0).
    enum class UForm
    {
        lich_fll = 1, // Lichtenstein (rotationally invariant) + FLL DC
        lich_amf = 2, // Lichtenstein (rotationally invariant) + AMF DC
        dud_fll  = 3, // Dudarev (simplified) + FLL DC -- default, implemented
    };

    Plus_U_Base();
    ~Plus_U_Base();

    /// allocate relevant data structures (base part, no LCAO types)
    void init_base(UnitCell& cell,
                   const int npol,
                   const int nspin,
                   const std::vector<int>& l_channel,
                   const bool yukawa_potential,
                   const double yukawa_lambda,
                   const std::string& global_readin_dir,
                   const std::string& global_out_dir,
                   const std::string& init_chg,
                   const std::string& device,
                   const std::vector<double>& hubbard_u,
                   const double uramping,
                   const int occ_mat_ctrl,
                   const int mixing_dftu);

    void uramping_update();
    bool u_converged();

    // u_current
    double get_u_current(int it) const { return u_current[it]; }
    void set_u_current(int it, double val) { u_current[it] = val; }
    int get_num_u_types() const { return static_cast<int>(u_current.size()); }

    // l_channel
    int get_l_channel(int it) const { return l_channel[it]; }
    bool has_l_channel(int it) const { return l_channel[it] != -1; }
    const std::vector<int>& get_l_channel_vec() const { return l_channel; }

    double get_uramping() const { return uramping; }
    int get_occ_mat_ctrl() const { return occ_mat_ctrl; }
    UForm get_form() const { return form; }


    /// Yukawa screening object (non-null only when use_yukawa())
    bool use_yukawa() const { return yukawa_ != nullptr; }
    YukawaScreening& yukawa() { return *yukawa_; }
    const YukawaScreening& yukawa() const { return *yukawa_; }

    // +U energy term
    double get_energy() const { return energy_u; }
    void set_energy(const double &e) { energy_u = e; }
    void set_double_energy() { energy_u *= 2.0; }


    /// get the U-term coefficient matrix pointer for the given spin channel
    /// (small (2l+1)^2 per-atom matrix, projector coefficient)
    ///
    /// nspin=1: isk is ignored, returns &uterm_mat[0]
    /// nspin=2: isk selects spin-up (0) or spin-down (1) half of the
    ///          split layout [all_up | all_dn]
    /// nspin=4: isk is ignored, returns &uterm_mat[0] (all Pauli blocks)
    const std::complex<double>* get_uterm_mat_spin(const int nspin, const int isk) const
    {
        if (nspin == 2 && isk == 1)
        {
            return uterm_mat.data() + uterm_mat.size() / 2;
        }
        return uterm_mat.data();
    }

    /// get size of the U-term coefficient matrix for a single spin channel
    ///
    /// nspin=1: full array size
    /// nspin=2: half of the total (one spin channel in split layout)
    /// nspin=4: full array size (all Pauli blocks are packed together)
    int get_size_uterm_mat_spin(const int nspin) const
    {
        return (nspin == 2) ? static_cast<int>(uterm_mat.size() / 2)
                            : static_cast<int>(uterm_mat.size());
    }

    int get_size_uterm_mat() const
    {
        return uterm_mat.size();
    }

    // dftu can be calculated only after occ_mat is ready
    bool is_occmat_ready() const { return occmat_ready_; }
    void set_occmat_ready() { occmat_ready_ = true; }
    void set_occmat_stale() { occmat_ready_ = false; }

    /// direct access to the occupation matrix object (new write path)
    OccupationMatrix& occmat() { return occmat_; }
    const OccupationMatrix& occmat() const { return occmat_; }

    /// access the occupation-matrix mixer (non-null only when mixing enabled)
    OccMatMixer& occ_mixer() { return *occ_mixer_; }
    const OccMatMixer& occ_mixer() const { return *occ_mixer_; }
    bool has_occ_mixer() const { return occ_mixer_ != nullptr; }

    // --- Accessors for free-function interfaces (e.g. DFTU_BASE::cal_occ_pw) ---
    const std::string& get_device() const { return device; }
    const std::vector<double>& get_u_current_vec() const { return u_current; }
    const std::vector<int>& get_uterm_mat_index() const { return uterm_mat_index; }
    std::vector<std::complex<double>>& get_uterm_mat() { return uterm_mat; }
    const std::vector<std::complex<double>>& get_uterm_mat() const { return uterm_mat; }
    double& energy_ref() { return energy_u; }

  private:
    // --- State flags ---
    // dftu can be calculated only after occ_mat is ready
    bool occmat_ready_ = false;

  protected:
    // --- U values and orbital configuration (set in init_base) ---
    std::vector<double> u_current;
    std::vector<double> u_target;
    std::vector<int> l_channel;

    double uramping = 0.0;
    int occ_mat_ctrl = 0;

    // --- Occupation matrices ---
    OccupationMatrix occmat_;

    // Occupation-matrix mixer; constructed only when mixing_dftu != 0.
    // Owns the flat uom/uom_save buffers and the mixing orchestration.
    std::unique_ptr<OccMatMixer> occ_mixer_;

    // --- Internal state ---
    double energy_u = 0.0;

    UForm form = UForm::dud_fll;
    std::string device;

    std::vector<std::complex<double>> uterm_mat;
    std::vector<int> uterm_mat_index;

    // Yukawa screening object; constructed only when use_yukawa() is true.
    // Owns the screening length, Slater integrals and derived U/J.
    std::unique_ptr<YukawaScreening> yukawa_;
};


#endif
