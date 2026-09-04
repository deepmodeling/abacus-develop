#ifndef OCC_MATRIX_H
#define OCC_MATRIX_H

#include "source_base/matrix.h"

#include <vector>

class UnitCell;

/**
 * @brief On-site occupation matrices for DFT+U.
 *
 * Owns the nested occ[iat][l][n][spin] matrices together with their saved
 * copy (used by mixing) and the iat->(l,n,m,ipol)->iwt lookup table.
 * Layout:
 *   nspin=1/2: occ[iat][l][n] has 2 spin channels of (2l+1)x(2l+1)
 *   nspin=4:   occ[iat][l][n] has 1 channel of (2l+1)*npol x (2l+1)*npol
 *              (all Pauli blocks packed together)
 */
class OccupationMatrix
{
  public:
    /// allocate occ/occ_save/iatlnmipol2iwt according to the cell
    void init(const UnitCell& cell,
              const std::vector<int>& orbital_corr,
              int nspin,
              int npol);

    // --- element access ---
    double get(int iat, int l, int n, int spin, int m1, int m2) const
    {
        return occ_[iat][l][n][spin](m1, m2);
    }
    double get_save(int iat, int l, int n, int spin, int m1, int m2) const
    {
        return occ_save_[iat][l][n][spin](m1, m2);
    }
    void set(int iat, int l, int n, int spin, int m1, int m2, double val)
    {
        occ_[iat][l][n][spin](m1, m2) = val;
    }

    /// direct matrix access for kernels that operate on whole blocks
    ModuleBase::matrix& mat(int iat, int l, int n, int spin)
    {
        return occ_[iat][l][n][spin];
    }
    const ModuleBase::matrix& mat(int iat, int l, int n, int spin) const
    {
        return occ_[iat][l][n][spin];
    }
    ModuleBase::matrix& mat_save(int iat, int l, int n, int spin)
    {
        return occ_save_[iat][l][n][spin];
    }
    const ModuleBase::matrix& mat_save(int iat, int l, int n, int spin) const
    {
        return occ_save_[iat][l][n][spin];
    }

    // --- bulk data access (used by IO and legacy call sites) ---
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& data() { return occ_; }
    const std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& data() const { return occ_; }
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& data_save() { return occ_save_; }
    const std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& data_save() const { return occ_save_; }

    // --- lookup table ---
    int iwt(int iat, int l, int n, int m, int ipol) const
    {
        return iatlnmipol2iwt_[iat][l][n][m][ipol];
    }
    const std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>& iatlnmipol2iwt() const
    {
        return iatlnmipol2iwt_;
    }

    // --- flat (de)serialization of one atom's correlated orbital ---
    /// nspin=1: fills occ with occ[iat][l][0][0] data
    /// nspin=2: fills occ with interleaved spin-up then spin-down data
    /// nspin=4: fills occ with occ[iat][l][0][0] data (all Pauli blocks)
    void get_flat(int iat, int l, std::vector<double>& occ) const;
    void set_flat(int iat, int l, int spin, const std::vector<double>& occ);

    // --- whole-array operations ---
    void zero(const UnitCell& cell, const std::vector<int>& orbital_corr);
    void copy_to_save(const UnitCell& cell, const std::vector<int>& orbital_corr);

    // --- flat mixing buffer (de)serialization over all atoms ---
    /// write occ into uom at offsets given by index (split spin layout)
    void write_to_flat(const UnitCell& cell,
                       const std::vector<int>& orbital_corr,
                       const std::vector<int>& index,
                       std::vector<double>& uom) const;
    /// read occ from uom at offsets given by index (split spin layout)
    void read_from_flat(const UnitCell& cell,
                        const std::vector<int>& orbital_corr,
                        const std::vector<int>& index,
                        const std::vector<double>& uom);
    /// write occ_save into uom_save (skips when uom_save is empty)
    void write_save_to_flat(const UnitCell& cell,
                            const std::vector<int>& orbital_corr,
                            const std::vector<int>& index,
                            std::vector<double>& uom_save) const;

    int nspin() const { return nspin_; }
    int npol() const { return npol_; }

  private:
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> occ_;
    std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>> occ_save_;
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> iatlnmipol2iwt_;
    int nspin_ = 0;
    int npol_ = 0;
};

namespace elecstate
{
/// occ = beta * occ + (1-beta) * occ_save on every atom's correlated orbital.
/// nspin-aware: nspin=4 mixes the single Pauli block, nspin=1/2 mixes both
/// spin channels. Replaces the duplicated LCAO k/gamma mixing loops.
void mix_occ_with_save(std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& occ_mat,
                       const std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>& occ_mat_save,
                       const UnitCell& cell,
                       const std::vector<int>& orbital_corr,
                       const int nspin,
                       const double beta);
} // namespace elecstate

#endif
