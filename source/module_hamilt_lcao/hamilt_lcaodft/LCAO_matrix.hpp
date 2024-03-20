#ifndef LCAO_MATRIX_HPP
#define LCAO_MATRIX_HPP

#include "LCAO_matrix.h"

//------------------------------------------------------
// DESCRIPTION:
// set 'dtype' matrix element (iw1_all, iw2_all) with 
// an input value 'v'
//------------------------------------------------------
template<typename T>
void LCAO_Matrix::set_mat2d(
    const int& global_ir, // index i for atomic orbital (row)
    const int& global_ic, // index j for atomic orbital (column)
    const T& v, // value for matrix element (i,j) 
    const Parallel_Orbitals& pv,
    T* HSloc)  //input pointer for store the matrix
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may be < 0 !!!
    const int ir = pv.global2local_row(global_ir);
    const int ic = pv.global2local_col(global_ic);

    const long index =
        ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()
        ? ic * static_cast<long>(pv.nrow) + ir
        : ir * static_cast<long>(pv.ncol) + ic;

    if (index >= pv.nloc)
    {
        std::cout << " iw1_all = " << global_ir << std::endl;
        std::cout << " iw2_all = " << global_ic << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " ParaV->nloc = " << pv.nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix", "set_mat2d");
    }	 

    //using input pointer HSloc
    HSloc[index] += v;
}

#endif