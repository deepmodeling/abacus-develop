#ifndef LCAO_HAMILT_H
#define LCAO_HAMILT_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Hamilt
{
    public:

    LCAO_Hamilt();
    ~LCAO_Hamilt();

    void grid_prepare(const Grid_Technique& gt, const ModulePW::PW_Basis& rhopw, const ModulePW::PW_Basis_Big& bigpw);

    // jingan add 2021-6-4
    void set_R_range_sparse();
    void calculate_STN_R_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_dSTN_R_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_STN_R_sparse_for_S(const double &sparse_threshold);
    void calculate_STN_R_sparse_for_T(const double &sparse_threshold);
    void calculat_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold);
    void calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold);
#ifdef __EXX
    template<typename Tdata> void calculate_HR_exx_sparse(
            const int &current_spin,
            const double &sparse_threshold,
            const int (&nmp)[3],
            const std::vector< std::map <int, std::map < std::pair<int, std::array<int,3>>, RI::Tensor<Tdata> > >>& Hexxs);
#endif
    void calculate_HSR_sparse(const int &current_spin, const double &sparse_threshold, const int (&nmp)[3]);
    void calculate_SR_sparse(const double &sparse_threshold);
    void clear_zero_elements(const int &current_spin, const double &sparse_threshold);
    void destroy_all_HSR_sparse(void);
    void calculate_TR_sparse(const double &sparse_threshold);
    void destroy_TR_sparse(void);
    void calculate_dH_sparse(const int &current_spin, const double &sparse_threshold);
    void destroy_dH_R_sparse(void);

    // used for gamma only algorithms.
    Gint_Gamma GG;

    // used for k-dependent grid integration.
    Gint_k GK;

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH genH;

    LCAO_Matrix* LM;
};

#endif
