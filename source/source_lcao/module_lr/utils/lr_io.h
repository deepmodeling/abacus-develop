#pragma once
#include "source_base/tool_title.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#ifdef __EXX
#include "source_lcao/module_ri/RI_Util.h" // for get_Born_von_Karmen_cells
#include <RI/global/Tensor.h>
#endif
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace LR_IO
{

inline void read_one_data(std::ifstream& ifs, double& data)
{
    std::string temp;
    ifs >> data >> temp;
}
inline void read_one_data(std::ifstream& ifs, std::complex<double>& data)
{
    double real, imag;
    ifs >> real >> imag;
    data = std::complex<double>(real, imag);
}
inline void read_one_data(double& re, double& im, std::complex<double>& out)
{
    out = std::complex<double>(re, im);
}
inline void read_one_data(double& re, double& im, double& out)
{
    out = re;
}
inline void set_zero_if_close(ModuleBase::Vector3<double>& vec, const double tol=1e-8)
{
    if (std::abs(vec.x) < tol) vec.x = 0.0;
    if (std::abs(vec.y) < tol) vec.y = 0.0;
    if (std::abs(vec.z) < tol) vec.z = 0.0;
}

void parse_band_out_file(const std::string& path, int& nbands_file, int& nk_file, int& nspin_file, int& nocc_file);

#ifdef __EXX
using TA = int;
using TC = std::array<int, 3>;
using TAC = std::pair<int, TC>;
template <typename T>
using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;

class RI_kRlist
{
  public:
    K_Vectors* klist = nullptr; // store fine kgrid if bse_use_fine_kgrid
    K_Vectors klist_coarse;
    TC period;
    std::vector<TC> Rlist;
    RI_kRlist() = default;
    RI_kRlist(const UnitCell& ucell, K_Vectors* pkv,
              const std::string& path, const int bse_use_fine_kgrid,
              const std::string& out_dir);
    ~RI_kRlist() = default;
    void read_kpts_coarse(const std::string& file, const UnitCell& ucell,
                          K_Vectors* const klist, const std::string& out_dir);
    void read_kpts_fine(const std::string& file, const UnitCell& ucell,
                        K_Vectors* const klist, const bool is_weighted,
                        const std::string& out_dir);
};

/// @brief vector as {ik, iband, <occ, ks_ene, gw_ene>}
/// @param ncore: as output, number of core orbitals parsed from file
std::vector<double> read_energy_qp(const int nocc,
                                   const int nvirt,
                                   const std::string& path,
                                   int& ncore,
                                   const int nk,
                                   const int nspin_tmp,
                                   const int nspin_file);

/// @brief same as read_energy_qp, but read from KS_band and GW_band files
std::vector<double> read_energy_qp_from_band_files(const K_Vectors& kv,
                                                   const int nocc,
                                                   const int nvirt,
                                                   int& ncore,
                                                   const std::string& path,
                                                   const int nk,
                                                   const int nspin_tmp,
                                                   const int nspin_file);

template <typename TK>
void read_librpa_eigenvectors(psi::Psi<TK>& wfc_ks,
                              psi::Psi<TK>& wfc_ks_global,
                              const std::string& path,
                              const int ncore,
                              const int nbands_file,
                              const int nspin_tmp,
                              const int nspin_file,
                              const int my_rank,
                              Parallel_Orbitals& pmat);

template <typename TK>
void read_librpa_eigenvectors_from_band_files(psi::Psi<TK>& wfc_ks,
                                              psi::Psi<TK>& wfc_ks_global,
                                              const std::string& path,
                                              const int ncore,
                                              const int nbands_file,
                                              const int nspin_tmp,
                                              const int nspin_file,
                                              const int my_rank,
                                              Parallel_Orbitals& pmat);

/// only for blocking by atom pairs (abacus type)
template <typename TCs, typename TR>
TLRI<TR> read_coulomb_mat_k(const std::string& path, const TLRI<TCs>& Cs, LR_IO::RI_kRlist& kRlist);

/// for any way of blocking (aims type)
template <typename TCs, typename TR>
TLRI<TR> read_coulomb_mat_general_k(const std::string& path, const TLRI<TCs>& Cs, LR_IO::RI_kRlist& kRlist);

/// @brief read Wxc(R) = Wc(R) + Vx(R) from file
template <typename Tdata, typename TR>
std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> read_Ws(const TLRI<TR>& Vs, const std::vector<TC>& Rlist);
#endif

} // namespace LR_IO
