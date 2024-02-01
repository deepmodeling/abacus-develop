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

// for debug
#ifdef _DEBUG
#include <unistd.h>

#include <cstring>
#include <fstream>

#include "src_pw/global.h"
#endif
// end debug

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
#ifdef _DEBUG
    char f_log[80];
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    std::ofstream log;
    if (myproc < 100)
    {
        sprintf(f_log, "transformer_%2.2d.log", myproc);
        log.open(f_log, std::ios::app);
        log << "start count nnz" << std::endl;
    }
    // count nonzeros value distribution of H and S
    static bool isCOUNTNONZERO = true;
    if (!isCOUNTNONZERO)
    {
        isCOUNTNONZERO = true;
        char plog_name[80];
        sprintf(plog_name, "HS_Distribution_%d.log", myproc);
        std::ofstream plog;
        plog.open(plog_name, std::ios::app);
        std::map<int, int> pH;
        countMatrixDistribution(nrow * ncol, H_2d, pH);
        std::map<int, int> pS;
        countMatrixDistribution(nrow * ncol, H_2d, pS);
        plog << "Element in H distribution:\n";
        // std::stringstream ss;
        // ss.str("");
        for (auto iter = pH.begin(); iter != pH.end(); ++iter)
        {
            // ss<<"p["<<iter->first<<"] : "<<iter->second<<std::endl;
            plog << "p[" << iter->first << "] : " << iter->second << std::endl;
        }
        // OUT(ofs_running,ss.str());
        // OUT(ofs_running, "Element in S distribution:");
        plog << "Element in S distribution:\n";
        // ss.str("");
        for (auto iter = pS.begin(); iter != pS.end(); ++iter)
        {
            // ss<<"p["<<iter->first<<"] : "<<iter->second<<std::endl;
            plog << "p[" << iter->first << "] : " << iter->second << std::endl;
        }
        // OUT(ofs_running,ss.str());
        plog.close();
    }
#endif

    int idx = 0;
    nnz = 0;
    colidx.clear();
    rowidx.clear();
#ifdef _DEBUG
    if (myproc < 100)
        log << "rowidx and colidx cleared" << std::endl;
#endif
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
#ifdef _DEBUG
        if (myproc < 100)
            log << "unknown layout: " << layout << std::endl;
#endif
        return 1;
    }
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "nnz is counted: " << nnz << std::endl;
        log.close();
    }
#endif
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
    // debug
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#ifdef _DEBUG
    std::ofstream log;
    if (myproc < 100)
    {
        char f_log[80];
        sprintf(f_log, "transformer_%2.2d.log", myproc);
        log.open(f_log, std::ios::app);
        log << "enter buildTransformParameter" << std::endl;
    }
#endif
    // end debug
    // count sender non-zeros elements
    sender_size = nnz;
    std::fill(sender_size_process.begin(), sender_size_process.end(), 0);
// debug
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "start translate ranks between group_data and group_trans" << std::endl;
        log << "sender_size (in BCD) = " << sender_size << std::endl;
    }
#endif
    // end debug
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

// debug
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "rank_data        rank_trans" << std::endl;
        for (int i = 0; i < nproc_data; ++i)
            log << i << "\t\t\t" << proc_map_data_trans[i] << std::endl;
    }
#endif
    // end debug

    for (int i = 0; i < nnz; ++i)
    {
        int l_col = colidx[i];
        int g_col = SRC_Matrix.globalCol(l_col);
        int dst_process;
        int dst_col = DST_Matrix.localCol(g_col, dst_process);
        int dst_process_trans = proc_map_data_trans[dst_process];
        /*
        // debug
        #ifdef _DEBUG
        log<<dst_process<<"\t\t";
        #endif
        // end debug
         MPI_Group_translate_ranks(DST_Matrix.group_data, 1, &dst_process,
                                   GROUP_TRANS, &dst_process_trans);
        // debug
        #ifdef _DEBUG
        log<<dst_process_trans<<std::endl;
        #endif
        // end debug
        */
        ++sender_size_process[dst_process_trans];
    }
// debug
#ifdef _DEBUG
    if (myproc < 100)
        log << "sender_size_process is creaated" << std::endl;
#endif
    // end debug

    // transfer sender index size to receiver index size
    MPI_Alltoall(&sender_size_process[0], 1, MPI_INT, &receiver_size_process[0], 1, MPI_INT, COMM_TRANS);
// debug
#ifdef _DEBUG
    if (myproc < 100)
        log << "receiver_size_process is got" << std::endl;
#endif
    // end debug

    // setup all2all parameters
    sender_displacement_process[0] = 0;
    for (int i = 1; i < NPROC_TRANS; ++i)
    {
        sender_displacement_process[i] = sender_displacement_process[i - 1] + sender_size_process[i - 1];
    }
// debug
#ifdef _DEBUG
    if (myproc < 100)
        log << "sender_displacement_process is creaated" << std::endl;
#endif
    // end debug

    receiver_displacement_process[0] = 0;
    receiver_size = receiver_size_process[0];
    for (int i = 1; i < NPROC_TRANS; ++i)
    {
        receiver_displacement_process[i] = receiver_displacement_process[i - 1] + receiver_size_process[i - 1];
        receiver_size += receiver_size_process[i];
    }
// debug
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "sender_size and receiver_displacement_process are creaated" << std::endl;
        log << "receiver_size (in CCS) = " << receiver_size << std::endl;
    }
#endif
    // end debug

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
// debug
#ifdef _DEBUG
    if (myproc < 100)
        log << "sender_index is got" << std::endl;
#endif
    // end debug

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
// debug
#ifdef _DEBUG
    if (myproc < 100)
        log << "receiver_index is got" << std::endl;
#endif
    // end debug

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
// debug
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "ccs parameter is built" << std::endl;
        log.close();
    }
#endif
    // end debug
    return 0;
}

int DistMatrixTransformer::newGroupCommTrans(DistBCDMatrix& SRC_Matrix,
                                             DistCCSMatrix& DST_Matrix,
                                             MPI_Group& GROUP_TRANS,
                                             MPI_Comm& COMM_TRANS)
{
// debug
#ifdef _DEBUG
    char f_log[80];
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    std::ofstream log;
    if (myproc < 100)
    {
        sprintf(f_log, "transformer_%2.2d.log", myproc);
        log.open(f_log, std::ios::app);
        // log<<std::endl<<"LOG of process: "<<myproc<<std::endl;
        log << "enter newGroupCommTrans" << std::endl;
    }
#endif
    // build transfortram communicator which contains both processes of BCD processors and
    // CCS processors with nonzero elements
    MPI_Group_union(DST_Matrix.get_group_data(), SRC_Matrix.get_group(), &GROUP_TRANS);
    MPI_Comm_create(MPI_COMM_WORLD, GROUP_TRANS, &COMM_TRANS);
// debug
#ifdef _DEBUG
    if (myproc < 100)
    {
        int trans_myid, trans_nproc;
        int trans_gid, trans_gproc;
        if (COMM_TRANS != MPI_COMM_NULL)
        {
            MPI_Comm_rank(COMM_TRANS, &trans_myid);
            MPI_Comm_size(COMM_TRANS, &trans_nproc);
        }
        else
        {
            trans_myid = -1;
            trans_nproc = -1;
            // trans_gid=-1;
            // trans_gproc=-1;
        }
        MPI_Group_rank(GROUP_TRANS, &trans_gid);
        MPI_Group_size(GROUP_TRANS, &trans_gproc);
        int BCD_myid, BCD_nproc;
        BCD_myid = SRC_Matrix.myproc;
        BCD_nproc = SRC_Matrix.nprocs;
        int BCD_gid, BCD_gproc;
        MPI_Group_rank(SRC_Matrix.group, &BCD_gid);
        MPI_Group_size(SRC_Matrix.group, &BCD_gproc);
        int CCS_myid, CCS_nproc;
        int CCS_gid, CCS_gproc;
        if (DST_Matrix.comm_data != MPI_COMM_NULL)
        {
            MPI_Comm_rank(DST_Matrix.comm_data, &CCS_myid);
            MPI_Comm_size(DST_Matrix.comm_data, &CCS_nproc);
        }
        else
        {
            CCS_myid = -1;
            CCS_nproc = -1;
            // CCS_gid=-1;
            // CCS_gproc=-1;
        }
        MPI_Group_rank(DST_Matrix.group_data, &CCS_gid);
        MPI_Group_size(DST_Matrix.group_data, &CCS_gproc);
        log << "myid in BCD:\t" << BCD_myid << "\tin CCS:\t" << CCS_myid << "\tin TRANS:\t" << trans_myid
            << "\tBCD_gid:\t" << BCD_gid << "\tCCS_gid:\t" << CCS_gid << "\ttrans_gid:\t" << trans_gid << std::endl;
        log << "nproc in BCD:\t" << BCD_nproc << "\tin CCS:\t" << CCS_nproc << "\tin TRANS:\t" << trans_nproc
            << "\tBCD_gproc:\t" << BCD_gproc << "\tCCS_gproc:\t" << CCS_gproc << "\ttrans_gproc:\t" << trans_gproc
            << std::endl;

        log << "COMM_TRANS is created" << std::endl;
        log.close();
    }
#endif
    // end debug
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
// debug
#ifdef _DEBUG
    char f_log[80];
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    std::ofstream log;
    if (myproc < 100)
    {
        sprintf(f_log, "transformer_%2.2d.log", myproc);
        log.open(f_log, std::ios::app);
        log << std::endl << "LOG of process: " << myproc << std::endl;
        log << "enter transformBCDtoCCS for H and S" << std::endl;
    }
#endif
    // end debug
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

#ifdef _DEBUG
        if (myproc < 100)
        {
            log << "nprocs: " << SRC_Matrix.nprocs << " ; myprow: " << SRC_Matrix.myprow
                << " ; mypcol: " << SRC_Matrix.mypcol << std::endl;
            log << "nblk:" << SRC_Matrix.nblk << " ; nrow: " << SRC_Matrix.get_nrow() << " ; ncol: " << SRC_Matrix.get_ncol()
                << std::endl;
            log << "layout:" << SRC_Matrix.get_layout() << std::endl;
            log << "ZERO = " << ZERO_Limit << std::endl;
            log << "DST_Matrix parameters:" << std::endl;
            log << "size: " << DST_Matrix.size << " ;nproc_data: " << DST_Matrix.nproc_data << std::endl;
            log << "start transforming H and S to CCS format" << std::endl;
        }
#endif
        // end debug

        // find out the non-zeros elements' positions
        std::vector<int> rowidx;
        std::vector<int> colidx;
        int nnz = 0;
#ifdef _DEBUG
        if (myproc < 100)
            log << "start counting nnz..." << std::endl;
#endif
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
#ifdef _DEBUG
        if (myproc < 100)
        {
            log << "NonZeroIndex is got, nnz is " << nnz << std::endl;
            log << "rowidx size: " << rowidx.size() << "; colidx size: " << colidx.size() << std::endl;
            /*
            if(SRC_Matrix.comm != MPI_COMM_NULL)
            {
                log<<"NonZeroIndex :"<<std::endl;
                if(SRC_Matrix.get_layout() == 'R' || SRC_Matrix.get_layout() == 'r')
                {
                    for(int i=0; i<nnz; ++i)
                    {
                        int HS_idx=rowidx[i]*SRC_Matrix.get_ncol()+colidx[i];
                        log<<rowidx[i]<<' '<<colidx[i]<<' '<<HS_idx;
                        log<<' '<<H_2d[HS_idx]<<' '<<S_2d[HS_idx]<<std::endl;
                    }
                }
                else
                {
                    for(int i=0; i<nnz; ++i)
                    {
                        int HS_idx=colidx[i]*SRC_Matrix.get_nrow()+rowidx[i];
                        log<<rowidx[i]<<' '<<colidx[i]<<' '<<HS_idx;
                        log<<' '<<H_2d[HS_idx]<<' '<<S_2d[HS_idx]<<std::endl;
                    }
                }
                log<<"nonzero index is output"<<std::endl;
            }
            else
            {
                log<<"no src_matrix elements in current process"<<std::endl;
            }
            */
        }
#endif

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
#ifdef _DEBUG
        if (myproc < 100)
            log << "Parameters are built" << std::endl;
#endif
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
#ifdef _DEBUG
        if (myproc < 100)
            log << "H sender_buffer is filled" << std::endl;
#endif
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
#ifdef _DEBUG
        if (myproc < 100)
            log << "H receiver_buffer is received" << std::endl;
#endif
        delete[] H_ccs;
        H_ccs = new double[receiver_size];
        buffer2CCSvalue(receiver_size, &buffer2ccsIndex[0], &receiver_buffer[0], H_ccs);
#ifdef _DEBUG
        if (myproc < 100)
            log << "H_ccs is received" << std::endl;
#endif

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
#ifdef _DEBUG
        if (myproc < 100)
            log << "S sender_buffer is filled" << std::endl;
#endif
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
#ifdef _DEBUG
        if (myproc < 100)
            log << "S receiver_buffer is received" << std::endl;
#endif
        delete[] S_ccs;
        S_ccs = new double[receiver_size];
        buffer2CCSvalue(receiver_size, &buffer2ccsIndex[0], &receiver_buffer[0], S_ccs);
#ifdef _DEBUG
        if (myproc < 100)
            log << "S_ccs is received" << std::endl;
#endif
    }
    // clear and return
    deleteGroupCommTrans(GROUP_TRANS, COMM_TRANS);
#ifdef _DEBUG
    if (myproc < 100)
    {
        log << "COMM_TRANS is deleted" << std::endl;
        log.close();
    }
#endif
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
// debug
#ifdef _DEBUG
    OUT(ofs_running, "transformCCStoBCD: start");
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // end debug
    int myproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
// debug
#ifdef _DEBUG
    std::ofstream log;
    if (myproc < 100)
    {
        char f_log[80];
        sprintf(f_log, "transformer_%2.2d.log", myproc);
        // MPI_Barrier(MPI_COMM_WORLD);
        log.open(f_log, std::ios::app);
        // MPI_Barrier(MPI_COMM_WORLD);
        log << "\nstart transform DMnzval to DM" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // end debug
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
#ifdef _DEBUG
        // MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "DM and EDM filled by 0" << std::endl;
// OUT(ofs_running, "transformCCStoBCD: DM and EDM filled by 0");
#endif
        // setup number of local elements to be transfered to each remote processes
        int NPROC_TRANS;
        MPI_Comm_size(COMM_TRANS, &NPROC_TRANS);
        // std::vector<int> sender_size_process(NPROC_TRANS);
        // std::vector<int> sender_displacement_process(NPROC_TRANS);
        // std::vector<int> receiver_size_process(NPROC_TRANS);
        // std::vector<int> receiver_displacement_process(NPROC_TRANS);
        int sender_size_process[NPROC_TRANS];
        int sender_displacement_process[NPROC_TRANS];
        int receiver_size_process[NPROC_TRANS];
        int receiver_displacement_process[NPROC_TRANS];
#ifdef _DEBUG
        if (myproc < 100)
            log << "NPROC_TRANS = " << NPROC_TRANS << std::endl;
        // MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "build process rank map from BCD to TRANS" << std::endl;
// OUT(ofs_running, "transformCCStoBCD: build process rank map from BCD to TRANS");
// MPI_Barrier(COMM_TRANS);
#endif
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

#ifdef _DEBUG
        // check process map from BCD comm to TRANS comm
        if (myproc < 100)
        {
            log << "check process map:\n";
            log << "pid in bcd\tpid in trans\n";
            for (int i = 0; i < nproc_bcd; ++i)
            {
                log << i << "\t\t" << proc_map_bcd_trans[i] << std::endl;
            }
            log << "check pid from prow and pcol int bcd to pid in trans\n";
            log << "p_row  p_col  p_bcd  p_trans\n";
            for (int i = 0; i < DST_Matrix.nprows; ++i)
            {
                for (int j = 0; j < DST_Matrix.npcols; ++j)
                {
                    int pid_bcd = DST_Matrix.pnum(i, j);
                    int pid_trans = proc_map_bcd_trans[pid_bcd];
                    log << i << "\t" << j << "\t" << pid_bcd << "\t" << pid_trans << std::endl;
                }
            }
            log << "setup alltoall parameters" << std::endl;
        }
        // OUT(ofs_running, "transformCCStoBCD: setup alltoall parameters");
        MPI_Barrier(COMM_TRANS);
#endif
        // setup sender_size_process
        // std::fill(sender_size_process.begin(), sender_size_process.end(), 0);
        for (int i = 0; i < NPROC_TRANS; ++i)
            sender_size_process[i] = 0;
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "sender_size_process is inited by 0" << std::endl;
        // OUT(ofs_running, "transformCCStoBCD: sender_size_process is inited by 0, size ", NPROC_TRANS);
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "display all columns and rows of nonzeros values:\n";
        int log_nnz = 0;
#endif
        for (int icol = 0; icol < SRC_Matrix.get_numcol_local(); ++icol)
        {
            int g_col = SRC_Matrix.globalCol(icol);
            int recv_pcol_bcd;
            int recv_col = DST_Matrix.localCol(g_col, recv_pcol_bcd);
            // #ifdef _DEBUG
            // log<<g_col<<"\n ";
            // #endif
            // OUT(ofs_running, "transformCCStoBCD: recv_pcol_bcd", recv_pcol_bcd);
            for (int rowidx = SRC_Matrix.get_colptr_local()[icol] - 1; rowidx < SRC_Matrix.get_colptr_local()[icol + 1] - 1; ++rowidx)
            {
                int g_row = SRC_Matrix.get_rowind_local()[rowidx] - 1;
                int recv_prow_bcd;
                int recv_row = DST_Matrix.localRow(g_row, recv_prow_bcd);
                int recv_proc_bcd = DST_Matrix.pnum(recv_prow_bcd, recv_pcol_bcd);
                int recv_proc = proc_map_bcd_trans[recv_proc_bcd];
                ++sender_size_process[recv_proc];
                // #ifdef _DEBUG
                // log<<" "<<g_row;
                // ++log_nnz;
                // #endif
            }
            // log<<"\n";
        }
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "sender_size_process is counted, total nonzeros are: " << log_nnz << std::endl;
            log << "target pid\tsize\n";
            for (int i = 0; i < NPROC_TRANS; i++)
            {
                log << i << "\t\t" << sender_size_process[i] << std::endl;
            }
        }
        // OUT(ofs_running, "transformCCStoBCD: sender_size_process is counted");
        MPI_Barrier(COMM_TRANS);
#endif

        // setup receiver_size_process
        // std::fill(receiver_size_process.begin(), receiver_size_process.end(), 0);
        for (int i = 0; i < NPROC_TRANS; ++i)
            receiver_size_process[i] = 0;
        MPI_Alltoall(&sender_size_process[0], 1, MPI_INT, &receiver_size_process[0], 1, MPI_INT, COMM_TRANS);
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "receiver_size_process is got" << std::endl;
            log << "target pid\tsize\n";
            for (int i = 0; i < NPROC_TRANS; i++)
            {
                log << i << "\t\t" << receiver_size_process[i] << std::endl;
            }
        }
// OUT(ofs_running, "transformCCStoBCD: receiver_size_process is got");
#endif

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
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "displacements are built" << std::endl;
            log << "check alltoallv parameters" << std::endl;
            for (int i = 0; i < NPROC_TRANS; ++i)
            {
                log << "pid_trans  sender_size_process  sender_displacement_process  receiver_size_process  "
                       "receiver_displacement_process"
                    << std::endl;
                log << i << "\t" << sender_size_process[i] << "\t\t\t" << sender_displacement_process[i] << "\t\t\t\t"
                    << receiver_size_process[i] << "\t\t\t" << receiver_displacement_process[i] << std::endl;
            }
        }
// OUT(ofs_running, "transformCCStoBCD: displacements are built");
#endif

        // setup up sender index and receiver index
        int sender_size = SRC_Matrix.get_nnzlocal();
        int* sender_index;
        double* sender_buffer;
        int* dst_index;
        int* receiver_index;
        double* receiver_buffer;
#ifdef _DEBUG
        if (myproc < 100)
        {
            log << "sender_size = " << sender_size << "; receiver_size = " << receiver_size << std::endl;
            log.flush();
            log << "start allocating sender_index, dst_index and receiver_index..." << std::endl;
            log.flush();
        }
#endif
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
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "; receiver_index size: ";
#endif
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
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "check BCD pnum" << std::endl;
            log.flush();
            for (int i = 0; i < DST_Matrix.nprows; ++i)
            {
                for (int j = 0; j < DST_Matrix.npcols; ++j)
                {
                    log << i << "\t" << j << "\t" << DST_Matrix.pnum(i, j) << std::endl;
                }
            }
            log << "source CCS matrix parameters:\n";
            log << "numColLocal: " << SRC_Matrix.numColLocal << std::endl;
            log << "pointer to beginning of each process is inited by sender_displacement_process" << std::endl;
            // log<<"icol"<<"\t"<<"g_col"<<"\t"<<"col(bcd)"<<"\t"<<"pcol(bcd)"<<std::endl;
            // log.flush();
        }
// MPI_Barrier(COMM_TRANS);
#endif

        int idx = 0;
#ifdef _DEBUG
        if (myproc < 100)
            log << "idx start at " << idx << std::endl;
#endif
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
#ifdef _DEBUG
                if (myproc < 100)
                {
                    if (recv_prow_bcd >= DST_Matrix.nprows || recv_prow_bcd < 0)
                    {
                        log << "ERROR: recv_prow_bcd error! recv_prow_bcd is " << recv_prow_bcd << "; max is "
                            << DST_Matrix.nprows << std::endl;
                        log.flush();
                    }
                }
#endif
                int recv_proc_bcd = DST_Matrix.pnum(recv_prow_bcd, recv_pcol_bcd);
#ifdef _DEBUG
                // MPI_Barrier(COMM_TRANS);
                if (myproc < 100)
                {
                    if (recv_proc_bcd > NPROC_TRANS || recv_proc_bcd < 0)
                    {
                        log << "ERROR: recv_proc_bcd outbound! recv_proc_bcd is " << recv_proc_bcd << "; max is "
                            << NPROC_TRANS << std::endl;
                        log.flush();
                    }
                }
#endif
                int recv_proc = proc_map_bcd_trans[recv_proc_bcd];
#ifdef _DEBUG
                // MPI_Barrier(COMM_TRANS);
                if (myproc < 100)
                {
                    if (p[recv_proc] >= sender_size || p[recv_proc] < 0)
                    {
                        log << "ERROR: sender_index's index outbound! " << std::endl;
                        log << recv_prow_bcd << " " << recv_pcol_bcd << recv_proc_bcd << " " << recv_proc << std::endl;
                        log << p[recv_proc] << " " << sender_size << std::endl;
                        log.flush();
                    }
                }
// MPI_Barrier(COMM_TRANS);
#endif
                sender_index[p[recv_proc]] = idx;
#ifdef _DEBUG
                // MPI_Barrier(COMM_TRANS);
                if (myproc < 100)
                {
                    if ((p[recv_proc] * 2 + 1) >= (2 * sender_size) || (p[recv_proc] * 2 + 1) < 0)
                    {
                        log << "ERROR: dst_index's index outbound! recv_proc:" << recv_proc
                            << "; p:" << p[recv_proc] * 2 + 1 << "; max is " << 2 * sender_size << std::endl;
                        log.flush();
                    }
                }
// MPI_Barrier(COMM_TRANS);
#endif
                dst_index[p[recv_proc] * 2] = recv_row;
                dst_index[p[recv_proc] * 2 + 1] = recv_col;
                ++p[recv_proc];
                ++idx;
            }
        }

#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        // check sender_index and dst_index
        if (myproc < 100)
        {
            for (int i = 0; i < sender_size; ++i)
            {
                if (sender_index[i] < 0 || sender_index[i] > SRC_Matrix.nnzLocal)
                {
                    log << "ERROR! sender_index outbound: " << i << " " << sender_index[i] << std::endl;
                    log.flush();
                }
            }
            for (int i = 0; i < 2 * sender_size; ++i)
            {
                if (dst_index[i] < 0 || dst_index[i] > DST_Matrix.size)
                {
                    log << "ERROR! dst_index outbound: " << i << " " << dst_index[i] << " " << DST_Matrix.size
                        << std::endl;
                    log.flush();
                }
            }
            log << "sender_index is built" << std::endl;
            log << "sender_size = " << sender_size << std::endl;
            // for(int i=0; i<sender_size; i+=sender_size/100)
            //     log<<i<<"\t"<<dst_index[2*i]<<"\t"<<dst_index[2*i+1]<<std::endl;
            // OUT(ofs_running, "transformCCStoBCD: sender_index is built");

            // save sender_index to file for debug
            /*std::ofstream log_sender_index;
            for(int i=0; i<NPROC_TRANS; ++i)
            {
                if(sender_size_process[i] > 0)
                {
                    sprintf(f_log, "sender_index_from_%2.2d_to_%2.2d.log", myproc_trans, i);
                    log_sender_index.open(f_log, std::ios::app);
                    for(int j=sender_displacement_process[i]; j<sender_displacement_process[i]+sender_size_process[i];
            ++j) log_sender_index<<sender_index[j]<<std::endl; log_sender_index.close();
                }
            }
            */
        }
#endif

        for (int i = 0; i < NPROC_TRANS; ++i)
        {
            sender_size_process[i] *= 2;
            sender_displacement_process[i] *= 2;
            receiver_size_process[i] *= 2;
            receiver_displacement_process[i] *= 2;
        }
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "Alltoall parameters for index array" << std::endl;
            log << "dst_index size:" << 2 * sender_size << "\t receiver_index size: " << 2 * receiver_size << std::endl;
            log << "pid_trans  sender_size_process  sender_displacement_process  receiver_size_process  "
                   "receiver_displacement_process"
                << std::endl;
            for (int i = 0; i < NPROC_TRANS; ++i)
            {
                log << i << "\t" << sender_size_process[i] << "\t\t" << sender_displacement_process[i] << "\t\t"
                    << receiver_size_process[i] << "\t\t" << receiver_displacement_process[i] << std::endl;
            }
            // save dst_index to file for debug
            /*std::ofstream log_dst_index;
            for(int i=0; i<NPROC_TRANS; ++i)
            {
                if(sender_size_process[i] > 0)
                {
                    sprintf(f_log, "dst_index_from_%2.2d_to_%2.2d.log", myproc_trans, i);
                    log_dst_index.open(f_log, std::ios::app);
                    for(int j=sender_displacement_process[i]; j<sender_displacement_process[i]+sender_size_process[i];
            ++j) log_dst_index<<dst_index[j]<<std::endl; log_dst_index.close();
                }
            }
            */
            log << "start alltoallv for index" << std::endl;
        }
        MPI_Barrier(COMM_TRANS);
// OUT(ofs_running, "transformCCStoBCD: sender_index is built");
#endif
        MPI_Alltoallv(&dst_index[0],
                      &sender_size_process[0],
                      &sender_displacement_process[0],
                      MPI_INT,
                      &receiver_index[0],
                      &receiver_size_process[0],
                      &receiver_displacement_process[0],
                      MPI_INT,
                      COMM_TRANS);
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
        {
            log << "receiver_index is got" << std::endl;
            log << "receiver_size is: " << receiver_size << std::endl;
            log.flush();
        }
/*
// save receiver_index to file for debug
std::ofstream log_rcv_index;
for(int i=0; i<NPROC_TRANS; ++i)
{
    log<<"receive index (from proc_trans "<<i<<") is from "<<receiver_displacement_process[i]<<" to
"<<receiver_displacement_process[i]+receiver_size_process[i]<<std::endl; if(receiver_size_process[i] > 0)
    {
        sprintf(f_log, "receiver_index_from_%2.2d_to_%2.2d.log", i, myproc_trans);
        log_rcv_index.open(f_log, std::ios::app);
        for(int j=receiver_displacement_process[i]; j<receiver_displacement_process[i]+receiver_size_process[i]; ++j)
            log_rcv_index<<receiver_index[j]<<std::endl;
        log_rcv_index.close();
    }
}
log<<"receiver_index values are saved"<<std::endl;
log.flush();
// MPI_Barrier(COMM_TRANS);

for(int i=0; i<receiver_size; ++i)
{
    if(receiver_index[i*2]<0)
    {
        log<<"ERROR! receiver_index(BCD)["<<2*i<<"] = "<<receiver_index[i*2]<<" < 0"<<std::endl;
        log.flush();
    }
    else if(receiver_index[i*2]>DST_Matrix.get_nrow())
    {
        log<<"ERROR! receiver_index(BCD)["<<2*i<<"] = "<<receiver_index[i*2]<<" > "<<DST_Matrix.get_nrow()<<std::endl;
        log.flush();
    }
    if(receiver_index[i*2+1]<0)
    {
        log<<"ERROR! receiver_index(BCD)["<<2*i+1<<"] = "<<receiver_index[i*2+1]<<" < 0"<<std::endl;
        log.flush();
    }
    else if(receiver_index[i*2+1]>DST_Matrix.get_ncol())
    {
        log<<"ERROR! receiver_index(BCD)["<<2*i+1<<"] = "<<receiver_index[i*2+1]<<" > "<<DST_Matrix.get_ncol()<<std::endl;
        log.flush();
    }
}
log<<"receiver_index values are checked"<<std::endl;
log.flush();
MPI_Barrier(COMM_TRANS);
// OUT(ofs_running, "transformCCStoBCD: receiver_index is got");
*/
#endif
        // reset size and displacement for transfering matrix value by alltoall
        for (int i = 0; i < NPROC_TRANS; ++i)
        {
            sender_size_process[i] /= 2;
            sender_displacement_process[i] /= 2;
            receiver_size_process[i] /= 2;
            receiver_displacement_process[i] /= 2;
        }
#ifdef _DEBUG
        if (myproc < 100)
        {
            log << "size_process and displacement_process are reset for buffer transform" << std::endl;
            log.flush();
        }
        MPI_Barrier(COMM_TRANS);
#endif

        // transfer DM
        // set up DM sender buffer
        for (int i = 0; i < sender_size; ++i)
        {
            sender_buffer[i] = DMnzvalLocal[sender_index[i]];
        }
#ifdef _DEBUG
        if (myproc < 100)
        {
            log << "DM(CCS) is put to sender_buffer" << std::endl;
            log.flush();
            // OUT(ofs_running, "transformCCStoBCD: DM(CCS) is put to sender_buffer");

            // check receiver_index, which may be changed after alltoall for buffer
            for (int i = 0; i < receiver_size; ++i)
            {
                if (receiver_index[i * 2] < 0)
                {
                    log << "ERROR! receiver_index(BCD)[" << 2 * i << "] = " << receiver_index[i * 2] << " < 0"
                        << std::endl;
                    log.flush();
                }
                else if (receiver_index[i * 2] > DST_Matrix.get_nrow())
                {
                    log << "ERROR! receiver_index(BCD)[" << 2 * i << "] = " << receiver_index[i * 2] << " > "
                        << DST_Matrix.get_nrow() << std::endl;
                    log.flush();
                }
                if (receiver_index[i * 2 + 1] < 0)
                {
                    log << "ERROR! receiver_index(BCD)[" << 2 * i + 1 << "] = " << receiver_index[i * 2 + 1] << " < 0"
                        << std::endl;
                    log.flush();
                }
                else if (receiver_index[i * 2 + 1] > DST_Matrix.get_ncol())
                {
                    log << "ERROR! receiver_index(BCD)[" << 2 * i + 1 << "] = " << receiver_index[i * 2 + 1] << " > "
                        << DST_Matrix.get_ncol() << std::endl;
                    log.flush();
                }
            }
            log << "receiver_index values are checked" << std::endl;
            log.flush();
            // check parameters for alltoall for buffer
            log << "pid_trans  sender_size_process  sender_displacement_process  receiver_size_process  "
                   "receiver_displacement_process"
                << std::endl;
            for (int i = 0; i < NPROC_TRANS; ++i)
            {
                log << i << "\t" << sender_size_process[i] << "\t\t" << sender_displacement_process[i] << "\t\t"
                    << receiver_size_process[i] << "\t\t" << receiver_displacement_process[i] << std::endl;
            }
            log.flush();
        }
        MPI_Barrier(COMM_TRANS);
#endif
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

#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "receiver_buffer is got from DM" << std::endl;
// OUT(ofs_running, "transformCCStoBCD: receiver_buffer is got from DM");
#endif
        // transform receiver_buffer to DM
        if (DST_Matrix.get_layout() == 'R' || DST_Matrix.get_layout() == 'r')
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = ix * DST_Matrix.get_ncol() + iy;
#ifdef _DEBUG
                if (myproc < 100)
                {
                    if (idx < 0 || idx >= DST_Matrix_elem)
                    {
                        log << "idx for DM ERROR: idx is " << idx << "; DM total size is " << DST_Matrix_elem
                            << std::endl;
                        log << "index number is " << 2 * i << " ix = " << ix << " iy = " << iy
                            << " ncol = " << DST_Matrix.get_ncol() << std::endl;
                        log.flush();
                    }
                }
#endif
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
#ifdef _DEBUG
                if (myproc < 100)
                {
                    if (idx < 0 || idx >= DST_Matrix_elem)
                    {
                        log << "idx for DM ERROR: idx is " << idx << "; DM total size is " << DST_Matrix_elem
                            << std::endl;
                        log << "index number is" << 2 * i << " ix = " << ix << " iy = " << iy
                            << " nrow = " << DST_Matrix.get_nrow() << std::endl;
                        log.flush();
                    }
                }
#endif
                DM[idx] = receiver_buffer[i];
            }
        }

#ifdef _DEBUG
        if (myproc < 100)
            log << "DM(BCD) is got from receiver_buffer" << std::endl;
        MPI_Barrier(COMM_TRANS);
// OUT(ofs_running, "transformCCStoBCD: DM(BCD) is got from receiver_buffer");
#endif
        // setup up sender buffer of EDM
        for (int i = 0; i < sender_size; ++i)
        {
            sender_buffer[i] = EDMnzvalLocal[sender_index[i]];
        }
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "EDM(CCS) is put to sender_buffer" << std::endl;
// OUT(ofs_running, "transformCCStoBCD: EDM(CCS) is put to sender_buffer");
#endif

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
#ifdef _DEBUG
        MPI_Barrier(COMM_TRANS);
        if (myproc < 100)
            log << "receiver_buffer is got from EDM" << std::endl;
// OUT(ofs_running, "transformCCStoBCD: receiver_buffer is got from EDM");
#endif
        // transform receiver_buffer to EDM
        if (DST_Matrix.get_layout() == 'R' || DST_Matrix.get_layout() == 'r')
        {
            int DST_Matrix_elem = DST_Matrix.get_nrow() * DST_Matrix.get_ncol();
            for (int i = 0; i < receiver_size; ++i)
            {
                int ix = receiver_index[2 * i];
                int iy = receiver_index[2 * i + 1];
                int idx = ix * DST_Matrix.get_ncol() + iy;
#ifdef _DEBUG
                if (myproc < 100)
                {
                    if (idx < 0 || idx >= DST_Matrix_elem)
                    {
                        log << "idx for EDM ERROR: idx is " << idx << "; EDM total size is " << DST_Matrix_elem
                            << std::endl;
                        log << "index number is" << 2 * i << " ix = " << ix << " iy = " << iy
                            << " ncol = " << DST_Matrix.get_ncol() << std::endl;
                        log.flush();
                    }
                }
#endif
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
#ifdef _DEBUG
                if (myproc < 100)
                {
                    if (idx < 0 || idx >= DST_Matrix_elem)
                    {
                        log << "idx for EDM ERROR: idx is " << idx << "; EDM total size is " << DST_Matrix_elem
                            << std::endl;
                        log << "index number is" << 2 * i << " ix = " << ix << " iy = " << iy
                            << " nrow = " << DST_Matrix.get_nrow() << std::endl;
                        log.flush();
                    }
                }
#endif
                EDM[idx] = receiver_buffer[i];
            }
        }
#ifdef _DEBUG
        if (myproc < 100)
            log << "EDM(BCD) is got from receiver_buffer" << std::endl;
        MPI_Barrier(COMM_TRANS);
#endif
        delete[] sender_index;
        delete[] sender_buffer;
        delete[] dst_index;
        delete[] receiver_index;
        delete[] receiver_buffer;
#ifdef _DEBUG
        if (myproc < 100)
            log << "work arrays are deleted" << std::endl;
#endif
    }
#ifdef _DEBUG
    if (myproc < 100)
        log << "OUT COMM_TRANS" << std::endl;
    if (myproc < 100)
        log << "before deleteGroupCommTrans" << std::endl;
#endif
    deleteGroupCommTrans(GROUP_TRANS, COMM_TRANS);
#ifdef _DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myproc < 100)
    {
        log << "COMM_TRANS is deleted" << std::endl;
        log.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    OUT(ofs_running, "transformCCStoBCD: finish job, COMM_TRANS is deleted");
#endif
    return 0;
}

} // namespace pexsi
#endif