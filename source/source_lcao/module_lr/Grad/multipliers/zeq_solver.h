#pragma once
#include "hamilt_zeq_left.h"
#include "hamilt_zeq_right.h"

namespace LR
{
    /// construct and solve Z-vector equation
    template<typename T>
    void Z_vector_equation(const T* const X,  ///< [in] transition amplitude X
        T* const Z,  ///< [out] lagrange multiplier Z
        const std::string& xc_kernel,
        const int& nstates,
        const int& nspin,
        const int& naos,
        const std::vector<int>& nocc,
        const std::vector<int>& nvirt,
        const UnitCell& ucell,
        const std::vector<double>& orb_cutoff,
        const Grid_Driver& gd,
        const psi::Psi<T>& psi_ks,
        const ModuleBase::matrix& eig_ks,
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri,
        const double& exx_alpha,
#endif 
        std::weak_ptr<PotHxcLR> pot,
        const K_Vectors& kv,
        const std::vector<Parallel_2D>& px,
        const Parallel_2D& pc,
        const Parallel_Orbitals& pmat,
        const std::string& spin_type = "singlet");
}

#include "zeq_solver.hpp"