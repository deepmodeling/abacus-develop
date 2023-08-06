#include "DM_gamma_2d_to_grid.h"
#include "module_base/global_function.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#ifdef __MPI
#include "module_base/blacs_connector.h"
#endif

DMgamma_2dtoGrid::DMgamma_2dtoGrid()
{
    sender_2D_index = nullptr;
    sender_size_process = nullptr;
    sender_displacement_process = nullptr;
    sender_buffer = nullptr;

    receiver_local_index = nullptr;
    receiver_size_process = nullptr;
    receiver_displacement_process = nullptr;
    receiver_buffer = nullptr;
}
DMgamma_2dtoGrid::~DMgamma_2dtoGrid()
{
    delete[] sender_2D_index;
    delete[] sender_size_process;
    delete[] sender_displacement_process;
    delete[] sender_buffer;

    delete[] receiver_local_index;
    delete[] receiver_size_process;
    delete[] receiver_displacement_process;
    delete[] receiver_buffer;
}
#ifdef __MPI
int DMgamma_2dtoGrid::setAlltoallvParameter(MPI_Comm comm_2D, int nbasis, int blacs_ctxt, int nblk, const int& loc_grid_dim, const int* global2local_grid)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "enter setAlltoallvParameter, nblk", nblk);
    ModuleBase::timer::tick("LOC", "Alltoall");
    this->comm_2D = comm_2D;
    // setup blacs parameters
    int nprows = 0;
    int npcols = 0;
    int nprocs = 0;
    int myprow = 0;
    int mypcol = 0;
    int myproc = 0;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

    Cblacs_pinfo(&myproc, &nprocs);
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nprocs",nprocs);

    size_t memory_sender = 0;
    size_t memory_receiver = 0;
    size_t memory_other = 0;
    // init data arrays
    delete[] sender_size_process;
    sender_size_process = new int[nprocs];
    delete[] sender_displacement_process;
    sender_displacement_process = new int[nprocs];
    //GlobalV::ofs_running << "checkpoint 2" << std::endl;
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lgd_now",lgd_now);
    memory_sender += nprocs * 2 * sizeof(int);

    receiver_size = loc_grid_dim * loc_grid_dim;
    delete[] receiver_size_process;
    receiver_size_process = new int[nprocs];
    delete[] receiver_displacement_process;
    receiver_displacement_process = new int[nprocs];
    delete[] receiver_local_index;
    receiver_local_index = new int[receiver_size];
    delete[] receiver_buffer;
    receiver_buffer = new double[receiver_size];

    memory_receiver += nprocs * 2 * sizeof(int) + receiver_size * (sizeof(int) + sizeof(double));

    int* trace_2D_row = new int[loc_grid_dim];
    int* trace_2D_col = new int[loc_grid_dim];
    int* trace_2D_prow = new int[loc_grid_dim];
    int* trace_2D_pcol = new int[loc_grid_dim];

    int* nRow_in_proc = new int[nprows];
    int* nCol_in_proc = new int[npcols];

    memory_other += (loc_grid_dim * 4 + nprows + npcols) * sizeof(int);

    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nprows",nprows);
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npcols",npcols);

    for (int i = 0; i < nprows; ++i)
    {
        nRow_in_proc[i] = 0;
    }
    for (int i = 0; i < npcols; ++i)
    {
        nCol_in_proc[i] = 0;
    }

    // count the number of elements to be received from each process
    auto localIndex = [](int globalindex, int nblk, int nprocs, int& myproc) -> int
        {
            myproc = int((globalindex % (nblk * nprocs)) / nblk);
            return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
        };
    for (int iGlobal = 0; iGlobal < nbasis; ++iGlobal)
    {
        int iLocalGrid = global2local_grid[iGlobal];
        if (iLocalGrid >= 0)
        {
            //trace_global[iLocalGrid]=iGlobal;
            int p;
            trace_2D_row[iLocalGrid] = localIndex(iGlobal, nblk, nprows, p);
            trace_2D_prow[iLocalGrid] = p;
            nRow_in_proc[trace_2D_prow[iLocalGrid]]++;
            trace_2D_col[iLocalGrid] = localIndex(iGlobal, nblk, npcols, p);
            trace_2D_pcol[iLocalGrid] = p;
            nCol_in_proc[trace_2D_pcol[iLocalGrid]]++;
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NLOCAL",nbasis);
    receiver_displacement_process[0] = 0;
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"receiver_displacement_process[0]",receiver_displacement_process[0]);
    for (int pnum = 0; pnum < nprocs; ++pnum)
    {
        int prow, pcol;
        Cblacs_pcoord(blacs_ctxt, pnum, &prow, &pcol);
        receiver_size_process[pnum] = nRow_in_proc[prow] * nCol_in_proc[pcol];

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pnum", pnum);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "prow", prow);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "pcol", pcol);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nRow_in_proc", nRow_in_proc[prow]);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nCol_in_proc", nCol_in_proc[pcol]);

        if (pnum > 0)
        {
            receiver_displacement_process[pnum] = receiver_displacement_process[pnum - 1] + receiver_size_process[pnum - 1];
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_size_process",receiver_size_process[nprocs-1]);

    // build the index to be received
    int* pos = new int[nprocs];
    int* receiver_2D_index = new int[receiver_size];
    for (int i = 0; i < nprocs; ++i)
    {
        pos[i] = receiver_displacement_process[i];
    }
    memory_other += nprocs * sizeof(int);
    memory_receiver += receiver_size * sizeof(int);

    for (int i = 0; i < loc_grid_dim; ++i)
    {
        int src_row = trace_2D_row[i];
        int src_prow = trace_2D_prow[i];
        for (int j = 0; j < loc_grid_dim; ++j)
        {
            int src_col = trace_2D_col[j];
            int src_idx = src_row * nbasis + src_col; // leanding dimension is set to nbasis for all processes

            int src_pcol = trace_2D_pcol[j];
            int src_proc = Cblacs_pnum(blacs_ctxt, src_prow, src_pcol);

            receiver_2D_index[pos[src_proc]] = src_idx;
            receiver_local_index[pos[src_proc]] = i * loc_grid_dim + j;
            ++pos[src_proc];
        }
    }
    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last receiver_2D_index",receiver_2D_index[lgd_now*lgd_now-1]);
    delete[] pos;
    delete[] trace_2D_row;
    delete[] trace_2D_col;
    delete[] trace_2D_prow;
    delete[] trace_2D_pcol;
    //delete[] trace_global;
    delete[] nRow_in_proc;
    delete[] nCol_in_proc;

    // send number of elements to be sent via MPI_Alltoall
    MPI_Alltoall(receiver_size_process, 1, MPI_INT,
        sender_size_process, 1, MPI_INT, comm_2D);

    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_size_process",sender_size_process[nprocs-1]);
    // setup sender buffer
    sender_size = sender_size_process[0];
    sender_displacement_process[0] = 0;
    for (int i = 1; i < nprocs; ++i)
    {
        sender_size += sender_size_process[i];
        sender_displacement_process[i] = sender_displacement_process[i - 1] + sender_size_process[i - 1];
    }

    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"sender_size",sender_size);
    delete[] sender_2D_index;
    sender_2D_index = new int[sender_size];
    delete[] sender_buffer;
    sender_buffer = new double[sender_size];
    memory_sender += sender_size * (sizeof(int) * sizeof(double));

    // send the index of the elements to be received via MPI_Alltoall
    MPI_Alltoallv(receiver_2D_index, receiver_size_process, receiver_displacement_process, MPI_INT,
        sender_2D_index, sender_size_process, sender_displacement_process, MPI_INT, comm_2D);


    GlobalV::ofs_running << "receiver_size is " << receiver_size << " ; receiver_size of each process is:\n";
    for (int i = 0; i < nprocs; ++i)
    {
        GlobalV::ofs_running << receiver_size_process[i] << " ";
    }
    GlobalV::ofs_running << std::endl;
    GlobalV::ofs_running << "sender_size is " << sender_size << " ; sender_size of each process is:\n";
    for (int i = 0; i < nprocs; ++i)
    {
        GlobalV::ofs_running << sender_size_process[i] << " ";
    }
    GlobalV::ofs_running << std::endl;

    // ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"last sender_2D_index",sender_2D_index[lgd_now*lgd_now-1]);
    delete[] receiver_2D_index;
    ModuleBase::Memory::record("LOC::A2A_receiv", memory_receiver);
    ModuleBase::Memory::record("LOC::A2A_sender", memory_sender);
    ModuleBase::Memory::record("LOC::A2A_other", memory_other);
    ModuleBase::timer::tick("LOC", "Alltoall");
    return 0;
}
#endif

// calculate the grid distributed DM matrix from 2D block-cyclic distributed DM matrix
// transform dm_gamma[is].c to this->DM[is]
void DMgamma_2dtoGrid::cal_dk_gamma_from_2D(
    std::vector<ModuleBase::matrix>& dm_gamma_2d,
    double*** dm_gamma_grid,
    int& nspin,
    int& nbasis,
    int& loc_grid_dim,
    std::ofstream& ofs_running)
{
    ModuleBase::timer::tick("LCAO_Charge", "dm_2dTOgrid");
#ifdef __DEBUG
    ModuleBase::GlobalFunc::OUT(ofs_running, "cal_dk_gamma_from_2D, NSPIN", NSPIN);
#endif

    for (int is = 0; is < nspin; ++is)
    {
        // put data from dm_gamma[is] to sender index
        int nNONZERO = 0;
        for (int i = 0; i < sender_size; ++i)
        {
            const int idx = sender_2D_index[i];
            const int icol = idx % nbasis;
            const int irow = (idx - icol) / nbasis;
            // sender_buffer[i]=wfc_dm_2d.dm_gamma[is](irow,icol);
            sender_buffer[i] = dm_gamma_2d[is](icol, irow); // sender_buffer is clomun major, 
            // so the row and column index should be switched
            if (sender_buffer[i] != 0) ++nNONZERO;
        }

#ifdef __DEBUG
        ModuleBase::GlobalFunc::OUT(ofs_running, "number of non-zero elements in sender_buffer", nNONZERO);
        ModuleBase::GlobalFunc::OUT(ofs_running, "sender_size", sender_size);
        ModuleBase::GlobalFunc::OUT(ofs_running, "last sender_buffer", sender_buffer[sender_size - 1]);
#endif

        // transform data via MPI_Alltoallv
#ifdef __MPI
        MPI_Alltoallv(sender_buffer, sender_size_process, sender_displacement_process, MPI_DOUBLE,
            receiver_buffer, receiver_size_process, receiver_displacement_process, MPI_DOUBLE, this->comm_2D);
#endif
        // put data from receiver buffer to this->DM[is]
        nNONZERO = 0;

        for (int i = 0; i < receiver_size; ++i)
        {
            const int idx = receiver_local_index[i];
            const int icol = idx % loc_grid_dim;
            const int irow = (idx - icol) / loc_grid_dim;
            dm_gamma_grid[is][irow][icol] = receiver_buffer[i];
            //DM[is][icol][irow]=receiver_buffer[i];
            if (receiver_buffer[i] != 0) ++nNONZERO;
        }

#ifdef __DEBUG
        ModuleBase::GlobalFunc::OUT(ofs_running, "number of non-zero elements in receiver_buffer", nNONZERO);
        ModuleBase::GlobalFunc::OUT(ofs_running, "receiver_size", receiver_size);
        ModuleBase::GlobalFunc::OUT(ofs_running, "last receiver_buffer", receiver_buffer[receiver_size - 1]);
#endif

    }
    ModuleBase::timer::tick("LCAO_Charge", "dm_2dTOgrid");
    return;
}