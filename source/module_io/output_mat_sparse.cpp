#include "output_mat_sparse.h"

#include "cal_r_overlap_R.h"
#include "write_HS_R.h"

namespace ModuleIO
{

    template<typename T>
    Output_Mat_Sparse<T>::Output_Mat_Sparse(
        int out_mat_hsR,
        int out_mat_dh,
        int out_mat_t,
        int out_mat_r,
        int istep,
        const ModuleBase::matrix& v_eff,
        const Parallel_Orbitals& pv,
        LCAO_Hamilt& UHM,
        Gint_k& gint_k, // mohan add 2024-04-01
        LCAO_Matrix& LM,
        const K_Vectors& kv,
        hamilt::Hamilt<T>* p_ham)
    : _out_mat_hsR(out_mat_hsR),
      _out_mat_dh(out_mat_dh),
      _out_mat_t(out_mat_t),
      _out_mat_r(out_mat_r),
      _istep(istep),
      _v_eff(v_eff),
      _pv(pv),
      _UHM(UHM),
      _gint_k(gint_k), // mohan add 2024-04-01
      _LM(LM),
      _kv(kv),
      _p_ham(p_ham)
{
}

template<>
void Output_Mat_Sparse<double>::write(void) 
{
}


template<>
void Output_Mat_Sparse<std::complex<double>>::write(void)
{
    if (_out_mat_hsR)
    {
        output_HS_R(_istep, this->_v_eff, this->_UHM, _kv, _p_ham);
    }

    if (_out_mat_t)
    {
        output_T_R(_istep, this->_UHM); // LiuXh add 2019-07-15
    }

    if (_out_mat_dh)
    {
		output_dH_R(
				_istep, 
				this->_v_eff, 
				this->_UHM, 
				this->_gint_k, // mohan add 2024-04-01
				this->_LM,
				_kv); // LiuXh add 2019-07-15
	}

    // add by jingan for out r_R matrix 2019.8.14
    if (_out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(this->_pv);
        if (_out_mat_hsR)
        {
            r_matrix.out_rR_other(_istep, this->_LM.output_R_coor);
        }
        else
        {
            r_matrix.out_rR(_istep);
        }
    }
}

template class Output_Mat_Sparse<double>;
template class Output_Mat_Sparse<std::complex<double>>;

} // namespace ModuleIO
