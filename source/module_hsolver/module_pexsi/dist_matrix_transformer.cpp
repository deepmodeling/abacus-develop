#ifdef __PEXSI
#include "dist_matrix_transformer.h"

#include <mpi.h>

#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "dist_bcd_matrix.h"
#include "dist_ccs_matrix.h"

namespace pexsi
{
// find the minimum index, the return value will be a non-negtive value index value if it is found, otherwise will be a
// negtive value the size_process and displacement_process array will be changed after the index is found isFirst:
// wether this function is called for the first time for a index array; nprocs: total number of processes size_process:
// the number of indices in each process displacement_process: the start position in each process index: the array
// contains the indices
inline int DistMatrixTransformer::MinimumIndexPosition(const bool isFirst,
                                                       const int nprocs,
                                                       int* size_process,
                                                       int* displacement_process,
                                                       const int* index)
{
    // usually the minimum index is continuous, so it will be a good idea to
    // check the one next to the previous index first.
    static int pre_position; // previous position in index array of minimum index,
    static int pre_process;  // the process contains previous index

    int minimum_index
        = INT_MAX; // the minimum index, initial value is a large number which is larger than any other index;
    int minimum_position = -1;
    int minimum_process = -1;

    if (isFirst)
    {
        for (int i = 0; i < nprocs; ++i)
        {
            if (size_process[i] > 0)
            {
                if (minimum_index > index[displacement_process[i]]) // find a smaller index
                {
                    minimum_position = displacement_process[i];
                    minimum_index = index[minimum_position];
                    minimum_process = i;
                }
            }
        }
        if (minimum_process >= 0) // find it!
        {
            ++displacement_process[minimum_process];
            --size_process[minimum_process];
        }
        pre_position = minimum_position;
        pre_process = minimum_process;
        return minimum_position;
    }
    else
    {
        // check the next one of pre_position
        if (size_process[pre_process] > 0 &&                    // the previous process still has elements
            index[pre_position + 1] == index[pre_position] + 1) // find it!
        {
            ++displacement_process[pre_process];
            --size_process[pre_process];
            ++pre_position;      // new pre_position is the next one
                                 // new pre_process keeps the same
            return pre_position; // current position is the new pre_position
        }

        // if the next one of pre_position is not the minimum one
        for (int i = 0; i < nprocs; ++i)
        {
            if (size_process[i] > 0)
            {
                if (minimum_index > index[displacement_process[i]])
                {
                    minimum_position = displacement_process[i];
                    minimum_index = index[minimum_position];
                    minimum_process = i;
                }
            }
        }
        if (minimum_process >= 0) // find it!
        {
            ++displacement_process[minimum_process];
            --size_process[minimum_process];
        }
        pre_position = minimum_position;
        pre_process = minimum_process;
        return minimum_position;
    }
}

inline void DistMatrixTransformer::buildCCSParameter(const int size,
                                                     const int nprocs,
                                                     std::vector<int> size_process,
                                                     std::vector<int> displacement_process,
                                                     const int* position_index,
                                                     DistCCSMatrix& DST_Matrix,
                                                     int* buffer2ccsIndex)
{
    // find the minimum one from left buffer index
    if (DST_Matrix.get_nnzlocal() <= 0)
        return;

    int pre_col = -1;
    int nnz_now = 0;
    int p_mini;
    p_mini = MinimumIndexPosition(true, nprocs, &size_process[0], &displacement_process[0], position_index);
    while (p_mini >= 0)
    {
        int index_mini = position_index[p_mini];
        int col_mini = index_mini / DST_Matrix.get_size(); //-DST_Matrix.firstCol;
        int row_mini = index_mini % DST_Matrix.get_size();
        if (col_mini > pre_col) // a new column starts, column pointer is a 1-based array
        {
            pre_col = col_mini;
            DST_Matrix.get_colptr_local()[col_mini] = nnz_now + 1;
        }
        DST_Matrix.get_rowind_local()[nnz_now] = row_mini + 1; // setup row index array, which is also 1-based
        // copy data from buffer to M, be careful M is a 0-based array
        buffer2ccsIndex[nnz_now] = p_mini;
        ++nnz_now;
        p_mini = MinimumIndexPosition(false, nprocs, &size_process[0], &displacement_process[0], position_index);
    }
    // The last element of colptrLocal is nnzLocal+1
    DST_Matrix.get_colptr_local()[DST_Matrix.get_numcol_local()] = nnz_now + 1;
}

inline void DistMatrixTransformer::buffer2CCSvalue(int nnzLocal,
                                                   int* buffer2ccsIndex,
                                                   double* buffer,
                                                   double* nzvalLocal)
{
    for (int i = 0; i < nnzLocal; ++i)
    {
        nzvalLocal[i] = buffer[buffer2ccsIndex[i]];
    }
}
inline void DistMatrixTransformer::countMatrixDistribution(int N, double* A, std::map<int, int>& P)
{
    for (int i = 0; i < N; ++i)
    {
        int key;
        if (fabs(A[i] < 1e-31))
            key = -100;
        else
            key = floor(log10(fabs(A[i])));
        ++P[key];
    }
}

// find out the index of non-zero elements
inline int DistMatrixTransformer::getNonZeroIndex(char layout,
                                                  const int nrow,
                                                  const int ncol,
                                                  double* H_2d,
                                                  double* S_2d,
                                                  const double ZERO_Limit,
                                                  int& nnz,
                                                  std::vector<int>& rowidx,
                                                  std::vector<int>& colidx)
{
    int idx = 0;
    nnz = 0;
    colidx.clear();
    rowidx.clear();
    if (layout == 'C' || layout == 'c')
    {
        for (int i = 0; i < ncol; ++i)
        {
            for (int j = 0; j < nrow; ++j)
            {
                idx = i * nrow + j;
                if (fabs(H_2d[idx]) > ZERO_Limit || fabs(S_2d[idx]) > ZERO_Limit)
                {
                    ++nnz;
                    colidx.push_back(i);
                    rowidx.push_back(j);
                }
            }
        }
    }
    else if (layout == 'R' || layout == 'r')
    {
        for (int i = 0; i < ncol; ++i)
        {
            for (int j = 0; j < nrow; ++j)
            {
                idx = j * ncol + i;
                if (fabs(H_2d[idx]) > ZERO_Limit || fabs(S_2d[idx]) > ZERO_Limit)
                {
                    ++nnz;
                    colidx.push_back(i);
                    rowidx.push_back(j);
                }
            }
        }
    }
    else
    {
        return 1;
    }
    return 0;
}

int DistMatrixTransformer::buildTransformParameter(DistBCDMatrix& SRC_Matrix,
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
                                                   std::vector<int>& buffer2ccsIndex)
{
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    sender_size = nnz;
    std::fill(sender_size_process.begin(), sender_size_process.end(), 0);
    // create process id map from group_data to group_trans
    int nproc_data;
    std::vector<int> proc_map_data_trans;
    if (myproc == 0)
    {
        MPI_Group_size(DST_Matrix.get_group_data(), &nproc_data);
        MPI_Bcast(&nproc_data, 1, MPI_INT, 0, COMM_TRANS);
        proc_map_data_trans.resize(nproc_data, 0);
        for (int i = 0; i < nproc_data; ++i)
        {
            MPI_Group_translate_ranks(DST_Matrix.get_group_data(), 1, &i, GROUP_TRANS, &proc_map_data_trans[i]);
        }
        MPI_Bcast(&proc_map_data_trans[0], nproc_data, MPI_INT, 0, COMM_TRANS);
    }
    else
    {
        MPI_Bcast(&nproc_data, 1, MPI_INT, 0, COMM_TRANS);
        proc_map_data_trans.resize(nproc_data, 0);
        MPI_Bcast(&proc_map_data_trans[0], nproc_data, MPI_INT, 0, COMM_TRANS);
    }

    for (int i = 0; i < nnz; ++i)
    {
        int l_col = colidx[i];
        int g_col = SRC_Matrix.globalCol(l_col);
        int dst_process;
        int dst_col = DST_Matrix.localCol(g_col, dst_process);
        int dst_process_trans = proc_map_data_trans[dst_process];
        ++sender_size_process[dst_process_trans];
    }
    // transfer sender index size to receiver index size
    MPI_Alltoall(&sender_size_process[0], 1, MPI_INT, &receiver_size_process[0], 1, MPI_INT, COMM_TRANS);
    // setup all2all parameters
    sender_displacement_process[0] = 0;
    for (int i = 1; i < NPROC_TRANS; ++i)
    {
        sender_displacement_process[i] = sender_displacement_process[i - 1] + sender_size_process[i - 1];
    }

    receiver_displacement_process[0] = 0;
    receiver_size = receiver_size_process[0];
    for (int i = 1; i < NPROC_TRANS; ++i)
    {
        receiver_displacement_process[i] = receiver_displacement_process[i - 1] + receiver_size_process[i - 1];
        receiver_size += receiver_size_process[i];
    }

    // setup receiver index
    // setup sender_index
    std::vector<int> sender_index(sender_size);
    for (int i = 0; i < nnz; ++i)
    {
        int l_col = colidx[i];
        int g_col = SRC_Matrix.globalCol(l_col);
        int dst_process;
        int dst_col = DST_Matrix.localCol(g_col, dst_process);
        int l_row = rowidx[i];
        int dst_row = SRC_Matrix.globalRow(l_row);
        sender_index[i] = dst_col * DST_Matrix.get_size() + dst_row;
    }

    // transfer index to receiver
    std::vector<int> receiver_index(receiver_size);
    MPI_Alltoallv(&sender_index[0],
                  &sender_size_process[0],
                  &sender_displacement_process[0],
                  MPI_INT,
                  &receiver_index[0],
                  &receiver_size_process[0],
                  &receiver_displacement_process[0],
                  MPI_INT,
                  COMM_TRANS);

    // setup buffer2ccsIndex based on receiver_index
    buffer2ccsIndex.resize(receiver_size);
    DST_Matrix.setnnz(receiver_size);
    buildCCSParameter(receiver_size,
                      NPROC_TRANS,
                      receiver_size_process,
                      receiver_displacement_process,
                      &receiver_index[0],
                      DST_Matrix,
                      &buffer2ccsIndex[0]);
    return 0;
}

int DistMatrixTransformer::newGroupCommTrans(DistBCDMatrix& SRC_Matrix,
                                             DistCCSMatrix& DST_Matrix,
                                             MPI_Group& GROUP_TRANS,
                                             MPI_Comm& COMM_TRANS)
{
    // build transfortram communicator which contains both processes of BCD processors and
    // CCS processors with nonzero elements
    MPI_Group_union(DST_Matrix.get_group_data(), SRC_Matrix.get_group(), &GROUP_TRANS);
    MPI_Comm_create(MPI_COMM_WORLD, GROUP_TRANS, &COMM_TRANS);
    return 0;
}

int DistMatrixTransformer::deleteGroupCommTrans(MPI_Group& GROUP_TRANS, MPI_Comm& COMM_TRANS)
{
    MPI_Group_free(&GROUP_TRANS);
    if (COMM_TRANS != MPI_COMM_NULL)
    {
        MPI_Comm_free(&COMM_TRANS);
    }
    return 0;
}

// transform two sparse matrices from block cyclic distribution (BCD) to Compressed Column Storage (CCS) distribution
// two destination matrices share the same non-zero elements positions
// if either of two elements in source matrices is non-zeros, the elements in the destination matrices are non-zero,
// even if one of them is acturely zero All matrices must have same MPI communicator
int DistMatrixTransformer::transformBCDtoCCS(DistBCDMatrix& SRC_Matrix,
                                             double* H_2d,
                                             double* S_2d,
                                             const double ZERO_Limit,
                                             DistCCSMatrix& DST_Matrix,
                                             double*& H_ccs,
                                             double*& S_ccs)
{
    MPI_Group GROUP_TRANS;
    MPI_Comm COMM_TRANS = MPI_COMM_NULL;
    newGroupCommTrans(SRC_Matrix, DST_Matrix, GROUP_TRANS, COMM_TRANS);
    if (COMM_TRANS != MPI_COMM_NULL)
    {
        // set up sender and receiver=
        int NPROC_TRANS;
        MPI_Comm_size(COMM_TRANS, &NPROC_TRANS);
        int sender_size;
        std::vector<int> sender_size_process(NPROC_TRANS);
        std::vector<int> sender_displacement_process(NPROC_TRANS);
        int receiver_size;
        std::vector<int> receiver_size_process(NPROC_TRANS);
        std::vector<int> receiver_displacement_process(NPROC_TRANS);

        // find out the non-zeros elements' positions
        std::vector<int> rowidx;
        std::vector<int> colidx;
        int nnz = 0;
        if (SRC_Matrix.get_comm() != MPI_COMM_NULL)
        {
            getNonZeroIndex(SRC_Matrix.get_layout(),
                            SRC_Matrix.get_nrow(),
                            SRC_Matrix.get_ncol(),
                            H_2d,
                            S_2d,
                            ZERO_Limit,
                            nnz,
                            rowidx,
                            colidx);
        }
        // build all2all transformation parameters and the map index of receiver buffer
        std::vector<int> buffer2ccsIndex;
        buildTransformParameter(SRC_Matrix,
                                DST_Matrix,
                                NPROC_TRANS,
                                GROUP_TRANS,
                                COMM_TRANS,
                                nnz,
                                rowidx,
                                colidx,
                                sender_size,
                                sender_size_process,
                                sender_displacement_process,
                                receiver_size,
                                receiver_size_process,
                                receiver_displacement_process,
                                buffer2ccsIndex);
// Do transformation
        std::vector<double> sender_buffer(sender_size);
        std::vector<double> receiver_buffer(receiver_size);
        // put H to sender buffer
        if (SRC_Matrix.get_layout() == 'R' || SRC_Matrix.get_layout() == 'r')
        {
            for (int i = 0; i < sender_size; ++i)
            {
                sender_buffer[i] = H_2d[rowidx[i] * SRC_Matrix.get_ncol() + colidx[i]];
            }
        }
        else
        {
            for (int i = 0; i < sender_size; ++i)
            {
                sender_buffer[i] = H_2d[colidx[i] * SRC_Matrix.get_nrow() + rowidx[i]];
            }
        }
        // do all2all transformation
        MPI_Alltoallv(&sender_buffer[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_DOUBLE,
                      &receiver_buffer[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_DOUBLE,
                      COMM_TRANS);
// collect H from receiver buffer
        delete[] H_ccs;
        H_ccs = new double[receiver_size];
        buffer2CCSvalue(receiver_size, &buffer2ccsIndex[0], &receiver_buffer[0], H_ccs);

        // put S to sender buffer
        if (SRC_Matrix.get_layout() == 'R' || SRC_Matrix.get_layout() == 'r')
        {
            for (int i = 0; i < sender_size; ++i)
            {
                sender_buffer[i] = S_2d[rowidx[i] * SRC_Matrix.get_ncol() + colidx[i]];
            }
        }
        else
        {
            for (int i = 0; i < sender_size; ++i)
            {
                sender_buffer[i] = S_2d[colidx[i] * SRC_Matrix.get_nrow() + rowidx[i]];
            }
        }
        // do all2all transformation
        MPI_Alltoallv(&sender_buffer[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_DOUBLE,
                      &receiver_buffer[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_DOUBLE,
                      COMM_TRANS);
// collect S from receiver buffer
        delete[] S_ccs;
        S_ccs = new double[receiver_size];
        buffer2CCSvalue(receiver_size, &buffer2ccsIndex[0], &receiver_buffer[0], S_ccs);
    }
    // clear and return
    deleteGroupCommTrans(GROUP_TRANS, COMM_TRANS);
    return 0;
}

// transform two sparse matrices from Compressed Column Storage (CCS) to block cyclic distribution (BCD) distribution
// two source matrices share the same non-zero elements positions
int DistMatrixTransformer::transformCCStoBCD(DistCCSMatrix& SRC_Matrix,
                                             double* DMnzvalLocal,
                                             double* EDMnzvalLocal,
                                             DistBCDMatrix& DST_Matrix,
                                             double* DM,
                                             double* EDM)
{
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    MPI_Group GROUP_TRANS;
    MPI_Comm COMM_TRANS = MPI_COMM_NULL;
    newGroupCommTrans(DST_Matrix, SRC_Matrix, GROUP_TRANS, COMM_TRANS);
    if (COMM_TRANS != MPI_COMM_NULL)
    {
        // init DM and EDM with 0
        for (int i = 0; i < DST_Matrix.get_nrow() * DST_Matrix.get_ncol(); ++i)
        {
            DM[i] = 0;
            EDM[i] = 0;
        }
        // setup number of local elements to be transfered to each remote processes
        int NPROC_TRANS;
        MPI_Comm_size(COMM_TRANS, &NPROC_TRANS);
        int sender_size_process[NPROC_TRANS];
        int sender_displacement_process[NPROC_TRANS];
        int receiver_size_process[NPROC_TRANS];
        int receiver_displacement_process[NPROC_TRANS];
        int nproc_bcd;
        std::vector<int> proc_map_bcd_trans;
        int myproc_trans;
        MPI_Comm_rank(COMM_TRANS, &myproc_trans);
        if (myproc_trans == 0)
        {
            MPI_Group_size(DST_Matrix.get_group(), &nproc_bcd);
            MPI_Bcast(&nproc_bcd, 1, MPI_INT, 0, COMM_TRANS);
            proc_map_bcd_trans.resize(nproc_bcd, 0);
            for (int i = 0; i < nproc_bcd; ++i)
            {
                MPI_Group_translate_ranks(DST_Matrix.get_group(), 1, &i, GROUP_TRANS, &proc_map_bcd_trans[i]);
            }
            MPI_Bcast(&proc_map_bcd_trans[0], nproc_bcd, MPI_INT, 0, COMM_TRANS);
        }
        else
        {
            MPI_Bcast(&nproc_bcd, 1, MPI_INT, 0, COMM_TRANS);
            proc_map_bcd_trans.resize(nproc_bcd, 0);
            MPI_Bcast(&proc_map_bcd_trans[0], nproc_bcd, MPI_INT, 0, COMM_TRANS);
        }
        // setup sender_size_process
        // std::fill(sender_size_process.begin(), sender_size_process.end(), 0);
        for (int i = 0; i < NPROC_TRANS; ++i)
            sender_size_process[i] = 0;
        for (int icol = 0; icol < SRC_Matrix.get_numcol_local(); ++icol)
        {
            int g_col = SRC_Matrix.globalCol(icol);
            int recv_pcol_bcd;
            int recv_col = DST_Matrix.localCol(g_col, recv_pcol_bcd);
    
            for (int rowidx = SRC_Matrix.get_colptr_local()[icol] - 1; rowidx < SRC_Matrix.get_colptr_local()[icol + 1] - 1; ++rowidx)
            {
                int g_row = SRC_Matrix.get_rowind_local()[rowidx] - 1;
                int recv_prow_bcd;
                int recv_row = DST_Matrix.localRow(g_row, recv_prow_bcd);
                int recv_proc_bcd = DST_Matrix.pnum(recv_prow_bcd, recv_pcol_bcd);
                int recv_proc = proc_map_bcd_trans[recv_proc_bcd];
                ++sender_size_process[recv_proc];
            }
            // log<<"\n";
        }
        // setup receiver_size_process
        // std::fill(receiver_size_process.begin(), receiver_size_process.end(), 0);
        for (int i = 0; i < NPROC_TRANS; ++i)
            receiver_size_process[i] = 0;
        MPI_Alltoall(&sender_size_process[0], 1, MPI_INT, &receiver_size_process[0], 1, MPI_INT, COMM_TRANS);
        // setup sender_displacement and receiver_displacement
        sender_displacement_process[0] = 0;
        receiver_displacement_process[0] = 0;
        int receiver_size = receiver_size_process[0];
        for (int i = 1; i < NPROC_TRANS; ++i)
        {
            sender_displacement_process[i] = sender_displacement_process[i - 1] + sender_size_process[i - 1];
            receiver_displacement_process[i] = receiver_displacement_process[i - 1] + receiver_size_process[i - 1];
            receiver_size += receiver_size_process[i];
        }

        // setup up sender index and receiver index
        int sender_size = SRC_Matrix.get_nnzlocal();
        int* sender_index;
        double* sender_buffer;
        int* dst_index;
        int* receiver_index;
        double* receiver_buffer;

        if (sender_size > 0)
        {
            sender_index = new int[sender_size];
            for (int i = 0; i < sender_size; ++i)
            {
                sender_index[i] = -1;
            }
            sender_buffer = new double[sender_size];
            dst_index = new int[2 * sender_size];
            for (int i = 0; i < 2 * sender_size; ++i)
            {
                dst_index[i] = -1;
            }
        }
        else
        {
            sender_index = new int[1];
            sender_index[0] = -1;
            sender_buffer = new double[1];
            dst_index = new int[2];
            dst_index[0] = -1;
            dst_index[1] = -1;
        }

        if (receiver_size > 0)
        {
            receiver_index = new int[2 * receiver_size];
            receiver_buffer = new double[receiver_size];
            for (int i = 0; i < 2 * receiver_size; ++i)
            {
                receiver_index[i] = -1;
            }
            for (int i = 0; i < receiver_size; ++i)
            {
                receiver_buffer[i] = -1;
            }
        }
        else
        {
            receiver_index = new int[2];
            receiver_buffer = new double[1];
            receiver_index[0] = -1;
            receiver_index[1] = -1;
            receiver_buffer[0] = -1;
        }

        // pointer to the first empty slot of each process
        // std::vector<int> p(sender_displacement_process);
        int p[NPROC_TRANS];
        for (int i = 0; i < NPROC_TRANS; ++i)
        {
            p[i] = sender_displacement_process[i];
        }

        int idx = 0;

        for (int icol = 0; icol < SRC_Matrix.get_numcol_local(); ++icol)
        {
            int g_col = SRC_Matrix.globalCol(icol);
            int recv_pcol_bcd;
            int recv_col = DST_Matrix.localCol(g_col, recv_pcol_bcd);
            for (int rowidx = SRC_Matrix.get_colptr_local()[icol] - 1; rowidx < SRC_Matrix.get_colptr_local()[icol + 1] - 1; ++rowidx)
            {
                int g_row = SRC_Matrix.get_rowind_local()[rowidx] - 1;
                int recv_prow_bcd;
                int recv_row = DST_Matrix.localRow(g_row, recv_prow_bcd);

                int recv_proc_bcd = DST_Matrix.pnum(recv_prow_bcd, recv_pcol_bcd);

                int recv_proc = proc_map_bcd_trans[recv_proc_bcd];

                sender_index[p[recv_proc]] = idx;

                dst_index[p[recv_proc] * 2] = recv_row;
                dst_index[p[recv_proc] * 2 + 1] = recv_col;
                ++p[recv_proc];
                ++idx;
            }
        }

        for (int i = 0; i < NPROC_TRANS; ++i)
        {
            sender_size_process[i] *= 2;
            sender_displacement_process[i] *= 2;
            receiver_size_process[i] *= 2;
            receiver_displacement_process[i] *= 2;
        }

        MPI_Alltoallv(&dst_index[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_INT,
                      &receiver_index[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_INT,
                      COMM_TRANS);

        // reset size and displacement for transfering matrix value by alltoall
        for (int i = 0; i < NPROC_TRANS; ++i)
        {
            sender_size_process[i] /= 2;
            sender_displacement_process[i] /= 2;
            receiver_size_process[i] /= 2;
            receiver_displacement_process[i] /= 2;
        }

        // transfer DM
        // set up DM sender buffer
        for (int i = 0; i < sender_size; ++i)
        {
            sender_buffer[i] = DMnzvalLocal[sender_index[i]];
        }

        // transfer sender buffer to receiver buffer
        MPI_Alltoallv(&sender_buffer[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_DOUBLE,
                      &receiver_buffer[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_DOUBLE,
                      COMM_TRANS);

        // transform receiver_buffer to DM
        if (DST_Matrix.get_layout() == 'R' || DST_Matrix.get_layout() == 'r')
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = ix * DST_Matrix.get_ncol() + iy;
                DM[idx] = receiver_buffer[i];
            }
        }
        else
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = iy * DST_Matrix.get_nrow() + ix;
                DM[idx] = receiver_buffer[i];
            }
        }

        // setup up sender buffer of EDM
        for (int i = 0; i < sender_size; ++i)
        {
            sender_buffer[i] = EDMnzvalLocal[sender_index[i]];
        }

        // transfer sender buffer to receiver buffer
        MPI_Alltoallv(&sender_buffer[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_DOUBLE,
                      &receiver_buffer[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_DOUBLE,
                      COMM_TRANS);

        // transform receiver_buffer to EDM
        if (DST_Matrix.get_layout() == 'R' || DST_Matrix.get_layout() == 'r')
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = ix * DST_Matrix.get_ncol() + iy;
                EDM[idx] = receiver_buffer[i];
            }
        }
        else
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = iy * DST_Matrix.get_nrow() + ix;
                EDM[idx] = receiver_buffer[i];
            }
        }

        delete[] sender_index;
        delete[] sender_buffer;
        delete[] dst_index;
        delete[] receiver_index;
        delete[] receiver_buffer;

    }
    deleteGroupCommTrans(GROUP_TRANS, COMM_TRANS);
    return 0;
}

} // namespace pexsi
#endif