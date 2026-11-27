#ifndef GET_WF_LCAO_H
#define GET_WF_LCAO_H

#include "source_base/parallel_grid.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

class Get_wf_lcao
{
  public:
    Get_wf_lcao(const psi::Psi<double>& psi, const Parallel_Orbitals& para_orb, int nspin);
    Get_wf_lcao(const psi::Psi<std::complex<double>>& psi, const Parallel_Orbitals& para_orb, int nspin);

    /// For gamma_only
    void begin_gamma(const UnitCell& ucell,
                     const Parallel_Grid& pgrid,
                     const std::vector<int>& out_wfc_norm,
                     const std::vector<int>& out_wfc_re_im,
                     const std::string& global_out_dir,
                     std::ofstream& ofs_running);

    /// For multi-k
    void begin_k(const UnitCell& ucell,
                 const Parallel_Grid& pgrid,
                 const K_Vectors& kv,
                 const std::vector<int>& out_wfc_norm,
                 const std::vector<int>& out_wfc_re_im,
                 const std::string& global_out_dir,
                 std::ofstream& ofs_running);

  private:
    const psi::Psi<double>* const psi_gamma_;
    const psi::Psi<std::complex<double>>* const psi_k_;
    const Parallel_Orbitals& para_orb_;
    const int nspin_;
    const int nbands_;

    void prepare_get_wf(std::ofstream& ofs_running);

    std::vector<int> select_bands(const std::vector<int>& out_wfc_kb) const;
};
#endif // GET_WF_LCAO_H
