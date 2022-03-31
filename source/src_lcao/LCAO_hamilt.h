#ifndef LCAO_HAMILT_H
#define LCAO_HAMILT_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "LCAO_gen_fixedH.h"
#include "gint_gamma.h"
#include "gint_k.h"

class LCAO_Hamilt
{
    public:

    LCAO_Hamilt();
    ~LCAO_Hamilt();

    void set_lcao_matrices(void);
        
    // used fro k-dependent Hamiltonian matrix.
    void calculate_Hk( const int &ik);
    
    // used for Gamma only Hamiltonian matrix.
    void calculate_Hgamma( const int &ik , vector<ModuleBase::matrix> dm_gamma);						// Peize Lin add ik 2016-12-03

    void calculate_STN_R(void); //LiuXh add 2019-07-15

    // jingan add 2021-6-4
    void set_R_range_sparse();
    void calculate_STN_R_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_STN_R_sparse_for_S(const double &sparse_threshold);
    void calculat_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold);
    void calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_HR_exx_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_HSR_sparse(const int &current_spin, const double &sparse_threshold);
    void calculate_SR_sparse(const double &sparse_threshold);
    void clear_zero_elements(const int &current_spin, const double &sparse_threshold);
    void destroy_all_HSR_sparse(void);

    // used for gamma only algorithms.
    Gint_Gamma GG;

    // used for k-dependent grid integration.
    Gint_k GK;

    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH genH;

    LCAO_Matrix* LM;
    
    // init S (overlap matrix) flag.
    bool init_s;

    private:

    // used for gamma only algorithms.
    void calculate_STNR_gamma(void);

    void calculate_STNR_gamma_B(void); //mohan add 2012-04-14

    void calculate_STNR_k(void);

};

#endif
