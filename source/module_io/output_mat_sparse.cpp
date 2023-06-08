#include "output_mat_sparse.h"

#include "cal_r_overlap_R.h"
#include "write_HS_R.h"

namespace ModuleIO
{

Output_Mat_Sparse::Output_Mat_Sparse(int out_mat_hsR,
                                     int out_mat_dh,
                                     int out_mat_t,
                                     int out_mat_r,
                                     bool is_md,
                                     int istep,
                                     int out_interval,
                                     const ModuleBase::matrix& v_eff,
                                     const Parallel_Orbitals& pv,
                                     LCAO_Hamilt& UHM,
                                     LCAO_Matrix& LM,
                                     const K_Vectors& kv)
    : _out_mat_hsR(out_mat_hsR),
      _out_mat_dh(out_mat_dh),
      _out_mat_t(out_mat_t),
      _out_mat_r(out_mat_r),
      _is_md(is_md),
      _istep(istep),
      _out_interval(out_interval),
      _v_eff(v_eff),
      _pv(pv),
      _UHM(UHM),
      _LM(LM),
      _kv(kv)
{
}

void Output_Mat_Sparse::write()
{
    if (_out_mat_hsR)
    {
        if (!_is_md || (_istep % _out_interval == 0))
        {
            output_HS_R(_istep, this->_v_eff, this->_UHM,
                        _kv); // LiuXh add 2019-07-15
        }                     // LiuXh add 2019-07-15
    }

    if (_out_mat_t)
    {
        if (!_is_md || (_istep % _out_interval == 0))
        {
            output_T_R(_istep, this->_UHM); // LiuXh add 2019-07-15
        }                                  // LiuXh add 2019-07-15
    }

    if (_out_mat_dh)
    {
        if (!_is_md || (_istep % _out_interval == 0))
        {
            output_dH_R(_istep, this->_v_eff, this->_UHM,
                        _kv); // LiuXh add 2019-07-15
        }                     // LiuXh add 2019-07-15
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

} // namespace ModuleIO