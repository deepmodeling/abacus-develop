#ifndef LCAO_DEEPKS_IO_H
#define LCAO_DEEPKS_IO_H

#ifdef __DEEPKS

#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_base/tool_title.h"

#include <string>
#include <torch/script.h>
#include <torch/torch.h>
#include <vector>

namespace LCAO_deepks_io
{

/// This file contains subroutines that contains interface with libnpy
/// since many arrays must be saved in numpy format
/// It also contains subroutines for printing density matrices
/// which is used in unit tests

/// There are 2 subroutines for printing and loading .npy file:
/// 1. print_dm : print density matrices
/// 2. load_npy_gedm : load gedm from .npy file

/// others print quantities in .npy format

/// 3. save_npy_d : descriptor -> deepks_dm_eig.npy
/// 4. save_npy_e : energy
/// 5. save_npy_f : force
/// 6. save_npy_gvx : gvx -> deepks_gradvx.npy
/// 7. save_npy_s : stress
/// 8. save_npy_gvepsl : gvepsl -> deepks_gvepsl.npy
/// 9. save_npy_o: orbital
/// 10. save_npy_orbital_precalc: orbital_precalc -> deepks_orbpre.npy
/// 11. save_npy_h : Hamiltonian
/// 12. save_npy_v_delta_precalc : v_delta_precalc -> deepks_vdpre.npy
/// 13. save_npy_phialpha : phialpha -> deepks_phialpha.npy
/// 14. save_npy_gevdm : grav_evdm -> deepks_gevdm.npy, can use phialpha and gevdm to calculate v_delta_precalc

/// print density matrices
template <typename TK>
void print_dm(const int nks, const int nlocal, const int nrow, const std::vector<std::vector<TK>>& dm);

void load_npy_gedm(const int nat, const int des_per_atom, double** gedm, double& e_delta, const int rank);

/// save descriptor
void save_npy_d(const int nat,
                const int des_per_atom,
                const int inlmax,
                const int* inl_l,
                const bool deepks_equiv,
                const std::vector<torch::Tensor>& d_tensor,
                const std::string& out_dir,
                const int rank);

// save energy
void save_npy_e(const double& e, /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/
                const std::string& e_file,
                const int rank);

// save force and gvx
void save_npy_f(const ModuleBase::matrix& f, /**<[in] \f$F_{base}\f$ or \f$F_{tot}\f$, in Ry/Bohr*/
                const std::string& f_file,
                const int nat,
                const int rank);

void save_npy_gvx(const int nat,
                  const int des_per_atom,
                  const torch::Tensor& gvx_tensor,
                  const std::string& out_dir,
                  const int rank);

// save stress and gvepsl
void save_npy_s(const ModuleBase::matrix& stress, /**<[in] \f$S_{base}\f$ or \f$S_{tot}\f$, in Ry/Bohr^3*/
                const std::string& s_file,
                const double& omega,
                const int rank);

void save_npy_gvepsl(const int nat,
                     const int des_per_atom,
                     const torch::Tensor& gvepsl_tensor,
                     const std::string& out_dir,
                     const int rank);

/// save orbital and orbital_precalc
void save_npy_o(const std::vector<double>& bandgap, /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/
                const std::string& o_file,
                const int nks,
                const int rank);

void save_npy_orbital_precalc(const int nat,
                              const int nks,
                              const int des_per_atom,
                              const torch::Tensor& orbital_precalc,
                              const std::string& out_dir,
                              const int rank);

// save Hamiltonian and v_delta_precalc(for deepks_v_delta==1)/phialpha+gevdm(for deepks_v_delta==2)
template <typename TK, typename TH>
void save_npy_h(const std::vector<TH>& hamilt,
                const std::string& h_file,
                const int nlocal,
                const int nks,
                const int rank);

template <typename TK>
void save_npy_v_delta_precalc(const int nat,
                              const int nks,
                              const int nlocal,
                              const int des_per_atom,
                              const torch::Tensor& v_delta_precalc_tensor,
                              const std::string& out_dir,
                              const int rank);

template <typename TK>
void save_npy_phialpha(const int nat,
                       const int nks,
                       const int nlocal,
                       const int inlmax,
                       const int lmaxd,
                       const torch::Tensor& phialpha_tensor,
                       const std::string& out_dir,
                       const int rank);

// Always real, no need for template now
void save_npy_gevdm(const int nat,
                    const int inlmax,
                    const int lmaxd,
                    const torch::Tensor& gevdm_tensor,
                    const std::string& out_dir,
                    const int rank);
}; // namespace LCAO_deepks_io

#endif
#endif
