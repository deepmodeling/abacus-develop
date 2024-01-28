#ifndef DISTMATRIXTRANSFORMER_H
#define DISTMATRIXTRANSFORMER_H

#include "dist_bcd_matrix.h"
#include "dist_ccs_matrix.h"
// transform a sparse matrix from block cyclic distribution (BCD) to Compressed Column Storage (CCS) distribution
// they should have same MPI communicator
// The local matrix of BCD is column-major order
// int transformBCDtoCCS(DistBCDMatrix &SRC_Matrix, double* H_2d, const double ZERO_Limit,
//                    DistCCSMatrix &DST_Matrix, double*& H_ccs);

// transform two sparse matrices from block cyclic distribution (BCD) to Compressed Column Storage (CCS) distribution
// two destination matrices share the same non-zero elements positions
// if either of two elements in source matrices is non-zeros, the elements in the destination matrices are non-zero,
// even if one of them is acturely zero All matrices must have same MPI communicator
namespace pexsi
{
int transformBCDtoCCS(DistBCDMatrix& SRC_Matrix,
                      double* H_2d,
                      double* S_2d,
                      const double ZERO_Limit,
                      DistCCSMatrix& DST_Matrix,
                      double*& H_ccs,
                      double*& S_ccs);

// int transformCCStoBCD(DistCCSMatrix& SRC_Matrix, double* DMnzvalLocal,
// DistBCDMatrix& DST_Matrix, double* DM_2d);

int transformCCStoBCD(DistCCSMatrix& SRC_Matrix,
                      double* DMnzvalLocal,
                      double* ENDnzvalLocal,
                      DistBCDMatrix& DST_Matrix,
                      double* DM_2d,
                      double* END_2d);
} // namespace pexsi
#endif // DISTMATRIXTRANSFORMER_H