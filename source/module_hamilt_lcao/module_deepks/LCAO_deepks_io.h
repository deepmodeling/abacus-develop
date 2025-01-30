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

/// 1. save_npy_d : descriptor -> deepks_dm_eig.npy
/// 2. save_npy_e : energy
/// 3. save_npy_h : Hamiltonian
/// 4. save_matrix2npy : ModuleBase::matrix -> .npy, for force, stress and orbital
/// 5. save_tensor2npy : torch::Tensor -> .npy, for precalculation variables

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
                const std::vector<torch::Tensor>& descriptor,
                const std::string& out_dir,
                const int rank);

// save energy
void save_npy_e(const double& e, /**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/
                const std::string& e_file,
                const int rank);

// save Hamiltonian
template <typename TK, typename TH>
void save_npy_h(const std::vector<TH>& hamilt,
                const std::string& h_file,
                const int nlocal,
                const int nks,
                const int rank);

void save_matrix2npy(const std::string& file_name,
                     const ModuleBase::matrix& matrix,
                     const int rank,
                     const double& scale = 1.0,
                     const char mode = 'N');

template <typename T>
void save_tensor2npy(const std::string& file_name, const torch::Tensor& tensor, const int rank);
}; // namespace LCAO_deepks_io

#endif
#endif
