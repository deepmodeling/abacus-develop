#include "output_mat_sparse.h"
#include "cal_r_overlap_R.h"
#include "write_HS_R.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h" // for ucell

namespace ModuleIO
{

    template<typename T>
    Output_Mat_Sparse<T>::Output_Mat_Sparse(
        int out_mat_hsR,
        int out_mat_dh,
        int out_mat_t,
        int out_mat_r,
        int istep,
        const ModuleBase::matrix &v_eff,
        const Parallel_Orbitals &pv,
        LCAO_gen_fixedH &gen_h, // mohan add 2024-04-02
        Gint_k &gint_k, // mohan add 2024-04-01
        LCAO_Matrix &lm,
        Grid_Driver &grid, // mohan add 2024-04-06
        const K_Vectors& kv,
        hamilt::Hamilt<T>* p_ham)
    : _out_mat_hsR(out_mat_hsR),
      _out_mat_dh(out_mat_dh),
      _out_mat_t(out_mat_t),
      _out_mat_r(out_mat_r),
      _istep(istep),
      _v_eff(v_eff),
      _pv(pv),
      _gen_h(gen_h), // mohan add 2024-04-02
      _gint_k(gint_k), // mohan add 2024-04-01
      _lm(lm),
      _grid(grid), // mohan add 2024-04-06
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
    //! generate a file containing the Hamiltonian and S(overlap) matrices 
    if (_out_mat_hsR)
    {
		output_HSR(
				_istep, 
				this->_v_eff, 
                this->_pv,
				this->_lm, 
                this->_grid,
				_kv, 
				_p_ham);
    }

    //! generate a file containing the kinetic energy matrix
    if (_out_mat_t)
    {
		output_TR(
				_istep, 
                GlobalC::ucell,
                this->_pv, 
				this->_lm, 
				this->_grid,
				this->_gen_h); // LiuXh add 2019-07-15
    }

    //! generate a file containing the derivatives of the Hamiltonian matrix (in Ry/Bohr)
    if (_out_mat_dh)
    {
		output_dHR(
				_istep, 
				this->_v_eff, 
                this->_gen_h,
				this->_gint_k, // mohan add 2024-04-01
				this->_lm,
                this->_grid, // mohan add 2024-04-06
				_kv); // LiuXh add 2019-07-15
	}

    // add by jingan for out r_R matrix 2019.8.14
    if (_out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(this->_pv);
        if (_out_mat_hsR)
        {
            r_matrix.out_rR_other(_istep, this->_lm.output_R_coor);
        }
        else
        {
            r_matrix.out_rR(_istep);
        }
    }

    return;
}

template class Output_Mat_Sparse<double>;
template class Output_Mat_Sparse<std::complex<double>>;

} // namespace ModuleIO
