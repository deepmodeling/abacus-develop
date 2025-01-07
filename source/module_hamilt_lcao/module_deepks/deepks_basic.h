#ifndef DEEPKS_BASIC_H
#define DEEPKS_BASIC_H

#ifdef __DEEPKS
#include "module_base/tool_title.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_basic.cpp
//------------------------

// The file contains 2 subroutines:
//  load_model : loads model for applying V_delta
//  cal_gevdm : d(des)/d(pdm), calculated using torch::autograd::grad

// load the trained neural network models
void load_model(const std::string& model_file, torch::jit::script::Module& model);

// calculate gevdm
void cal_gevdm(const int nat,
               const int inlmax,
               const int* inl_l,
               const std::vector<torch::Tensor>& pdm,
               std::vector<torch::Tensor>& gevdm);

} // namespace DeePKS_domain
#endif
#endif