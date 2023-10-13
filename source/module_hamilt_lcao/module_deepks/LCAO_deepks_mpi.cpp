//This file contains only one subroutine, allsum_deepks
//which is used to perform allsum on a two-level pointer
//It is used in a few places in the deepks code

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"

#ifdef __MPI
void LCAO_Deepks::allsum_deepks(
    int inlmax, //first dimension
    int ndim, //second dimension
    double** mat) //the array being reduced 
{
    for(int inl=0;inl<inlmax;inl++)
    {
        Parallel_Reduce::reduce_all(mat[inl], ndim);
    }
}
#endif

#endif