#ifndef DEEPKS_BASIC_H
#define DEEPKS_BASIC_H

#ifdef __DEEPKS
#include "LCAO_deepks_io.h"
#include "module_base/parallel_reduce.h"
#include "module_base/tool_title.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_basic.cpp
//------------------------

// The file contains 2 subroutines:
// 1. load_model : loads model for applying V_delta
// 2. cal_gevdm : d(des)/d(pdm), calculated using torch::autograd::grad
// 3. cal_edelta_gedm : calculates E_delta and d(E_delta)/d(pdm)
//       this is the term V(D) that enters the expression H_V_delta = |alpha>V(D)<alpha|
//       caculated using torch::autograd::grad
// 4. check_gedm : prints gedm for checking
// 5. cal_edelta_gedm_equiv : calculates E_delta and d(E_delta)/d(pdm) for equivariant version

// load the trained neural network models
void load_model(const std::string& model_file, torch::jit::script::Module& model);

// calculate gevdm
void cal_gevdm(const int nat,
               const int inlmax,
               const int* inl_l,
               const std::vector<torch::Tensor>& pdm,
               std::vector<torch::Tensor>& gevdm);

/// calculate partial of energy correction to descriptors
void cal_edelta_gedm(const int nat,
              const int lmaxd,
              const int nmaxd,
              const int inlmax,
              const int des_per_atom,
              const int* inl_l,
              const std::vector<torch::Tensor>& descriptor,
              const std::vector<torch::Tensor>& pdm,
              torch::jit::script::Module& model_deepks,
              double** gedm,
              double& E_delta);
void check_gedm(const int inlmax, const int* inl_l, double** gedm);
void cal_edelta_gedm_equiv(const int nat,
                    const int lmaxd,
                    const int nmaxd,
                    const int inlmax,
                    const int des_per_atom,
                    const int* inl_l,
                    const std::vector<torch::Tensor>& descriptor,
                    double** gedm,
                    double& E_delta);

} // namespace DeePKS_domain
#endif
#endif