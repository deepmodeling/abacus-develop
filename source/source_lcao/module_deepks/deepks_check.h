#ifndef DEEPKS_CHECK_H
#define DEEPKS_CHECK_H

#ifdef __MLALGO

#include <string>
#include <torch/script.h>
#include <torch/torch.h>

namespace DeePKS_domain
{
//------------------------
// deepks_check.cpp
//------------------------

// This file contains subroutines for checking files

// There are 1 subroutines in this file:
// 1. check_tensor, which is used for tensor data checking

template <typename T>
void check_tensor(const torch::Tensor& tensor, const std::string& filename, const int rank);

} // namespace DeePKS_domain

#endif
#endif
