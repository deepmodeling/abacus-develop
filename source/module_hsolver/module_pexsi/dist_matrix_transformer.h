#ifndef DISTMATRIXTRANSFORMER_H
#define DISTMATRIXTRANSFORMER_H

#include <mpi.h>

#include <vector>
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
class DistBCDMatrix;
class DistCCSMatrix;

namespace DistMatrixTransformer
{
int MinimumIndexPosition(const bool isFirst,
                         const int nprocs,
                         int* size_process,
                         int* displacement_process,
                         const int* index);

void buildCCSParameter(const int size,
                       const int nprocs,
                       std::vector<int> size_process,
                       std::vector<int> displacement_process,
                       const int* position_index,
                       DistCCSMatrix& DST_Matrix,
                       int* buffer2ccsIndex);

void buffer2CCSvalue(int nnzLocal, int* buffer2ccsIndex, double* buffer, double* nzvalLocal);

void countMatrixDistribution(int N, double* A, std::map<int, int>& P);

int getNonZeroIndex(char layout,
                    const int nrow,
                    const int ncol,
                    double* H_2d,
                    double* S_2d,
                    const double ZERO_Limit,
                    int& nnz,
                    std::vector<int>& rowidx,
                    std::vector<int>& colidx);

int buildTransformParameter(DistBCDMatrix& SRC_Matrix,
                            DistCCSMatrix& DST_Matrix,
                            const int NPROC_TRANS,
                            MPI_Group& GROUP_TRANS,
                            MPI_Comm& COMM_TRANS,
                            const int nnz,
                            std::vector<int>& rowidx,
                            std::vector<int>& colidx,
                            int& sender_size,
                            std::vector<int>& sender_size_process,
                            std::vector<int>& sender_displacement_process,
                            int& receiver_size,
                            std::vector<int>& receiver_size_process,
                            std::vector<int>& receiver_displacement_process,
                            std::vector<int>& buffer2ccsIndex);

int newGroupCommTrans(DistBCDMatrix& SRC_Matrix,
                      DistCCSMatrix& DST_Matrix,
                      MPI_Group& GROUP_TRANS,
                      MPI_Comm& COMM_TRANS);

int deleteGroupCommTrans(MPI_Group& GROUP_TRANS, MPI_Comm& COMM_TRANS);

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
                      double* ED_2d);
}; // namespace DistMatrixTransformer
} // namespace pexsi
#endif // DISTMATRIXTRANSFORMER_H