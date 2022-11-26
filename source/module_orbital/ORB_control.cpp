#include "ORB_control.h"
#include "ORB_gen_tables.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_common.h"
#include "../src_io/wf_local.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/memory.h"

//#include "build_st_pw.h"

ORB_control::ORB_control(
    const bool& gamma_only_in,
    const int& nlocal_in,
    const int& nbands_in,
    const int& nspin_in,
    const int& dsize_in,
    const int& nb2d_in,
    const int& dcolor_in,
    const int& drank_in,
    const int& myrank_in,
    const std::string& calculation_in,
    const std::string& ks_solver_in) :
    gamma_only(gamma_only_in),
    nlocal(nlocal_in),
    nbands(nbands_in),
    nspin(nspin_in),
    dsize(dsize_in),
    nb2d(nb2d_in),
    dcolor(dcolor_in),
    drank(drank_in),
    myrank(myrank_in),
    calculation(calculation_in),
    ks_solver(ks_solver_in),
    setup_2d(true)
{
    this->ParaV.nspin = nspin_in;
}

ORB_control::ORB_control() :
    setup_2d(false)
{}
ORB_control::~ORB_control()
{}

void ORB_control::read_orb_first(
    std::ofstream& ofs_in,
    LCAO_Orbitals& orb,
    const int& ntype, // mohan add 2021-04-26
    const int& lmax, // mohan add 2021-04-26 
    const double& lcao_ecut_in, // mohan add 2021-04-16
    const double& lcao_dk_in, // mohan add 2021-04-16
    const double& lcao_dr_in, // mohan add 2021-04-16
    const double& lcao_rmax_in, // mohan add 2021-04-16
    const bool& deepks_setorb,
    const int& out_mat_r,
    const bool& force_flag, // mohan add 2021-05-07
    const int& my_rank // mohan add 2021-04-26
)
{
    ModuleBase::TITLE("ORB_control", "read_orb_first");
    ModuleBase::timer::tick("ORB_control", "read_orb_first");

    /////////////////////////////////////////////////////////////////
    /// (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
    ///
    /// (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
    /////////////////////////////////////////////////////////////////

    // mohan add 2021-04-16
    assert(ntype > 0);
    assert(lmax >= 0);
    assert(lcao_ecut_in > 0.0);
    assert(lcao_dk_in > 0.0);
    assert(lcao_dr_in > 0.0);
    assert(lcao_rmax_in > 0.0);

    // mohan add 2021-04-16
    orb.ecutwfc = lcao_ecut_in;
    orb.dk = lcao_dk_in;
    orb.dR = lcao_dr_in;
    orb.Rmax = lcao_rmax_in;

    orb.Read_Orbitals(
        ofs_in,
        ntype,
        lmax,
        deepks_setorb,
        out_mat_r,
        force_flag,
        my_rank);

    ModuleBase::timer::tick("ORB_control", "read_orb_first");
    return;
}

void ORB_control::set_orb_tables(
    std::ofstream& ofs_in,
    ORB_gen_tables& OGT,
    LCAO_Orbitals& orb,
    const double& lat0,
    const bool& deepks_setorb,
    const int& Lmax_exx,
    const int& nprojmax,
    const int* nproj,
    const Numerical_Nonlocal* beta_)
{
    ModuleBase::TITLE("ORB_control", "set_orb_tables");
    ModuleBase::timer::tick("ORB_control", "set_orb_tables");


#ifdef __NORMAL

#else
    if (calculation == "test")
    {
        ModuleBase::timer::tick("ORB_control", "set_orb_tables");
        return;
    }
#endif


    ///////////////////////////////////////////////////////////////////
    /// (2) FUNCTION : Generate Gaunt_Coefficients and S-table using OGT.init
    /// 	   Must have 'Numerical Orbital' infomation
    ///
    /// (2) RESULT : we have tabulated S table for use.
    ///////////////////////////////////////////////////////////////////
    const int job0 = 3;
    /// job0 :
    /// 1. generate overlap table
    /// 2. generate kinetic table
    /// 3. generate overlap & kinetic table
    OGT.gen_tables(ofs_in, job0, orb, Lmax_exx, deepks_setorb, nprojmax, nproj, beta_);
    // init lat0, in order to interpolated value from this table.

    assert(lat0 > 0.0);
    OGT.set_unit(lat0);

    ModuleBase::timer::tick("ORB_control", "set_orb_tables");
    return;
}

void ORB_control::clear_after_ions(
    ORB_gen_tables& OGT,
    LCAO_Orbitals& orb,
    const bool& deepks_setorb,
    const int* nproj_)
{
    ModuleBase::TITLE("ORB_control", "clear_after_ions");
    OGT.MOT.Destroy_Table(orb);
    OGT.tbeta.Destroy_Table_Beta(orb.get_ntype(), orb.Phi, nproj_);

    //caoyu add 2021-03-18
    if (deepks_setorb)
    {
        OGT.talpha.Destroy_Table_Alpha(orb);
    }
    return;
}


void ORB_control::setup_2d_division(std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "setup_2d_division");
    ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << std::endl;

    // (1) calculate nrow, ncol, nloc.
    if (ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver")
    {
        ofs_running << " divide the H&S matrix using 2D block algorithms." << std::endl;
#ifdef __MPI
        // storage form of H and S matrices on each processor
        // is determined in 'divide_HS_2d' subroutine
        this->divide_HS_2d(DIAG_WORLD, ofs_running, ofs_warning);
#else
        ModuleBase::WARNING_QUIT("LCAO_Matrix::init", "diago method is not ready.");
#endif
    }
    else
    {
        // the full matrix
        this->ParaV.nloc = nlocal * nlocal;
    }

    // (2) set the trace, then we can calculate the nnr.
    // for 2d: calculate po.nloc first, then trace_loc_row and trace_loc_col
    // for O(N): calculate the three together.
    this->set_trace(ofs_running);
}


void ORB_control::set_parameters(std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "set_parameters");

    Parallel_Orbitals* pv = &this->ParaV;
    // set loc_size
    if (gamma_only)//xiaohui add 2014-12-21
    {
        pv->loc_size = nbands / dsize;

        // mohan add 2012-03-29
        if (pv->loc_size == 0)
        {
            ofs_warning << " loc_size=0" << " in proc " << myrank + 1 << std::endl;
            ModuleBase::WARNING_QUIT("ORB_control::set_parameters", "nbands < ncpus");
        }

        if (drank < nbands % dsize) pv->loc_size += 1;
        if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "local size", pv->loc_size);

        // set loc_sizes
        delete[] pv->loc_sizes;
        pv->loc_sizes = new int[dsize];
        ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, dsize);

        pv->lastband_in_proc = 0;
        pv->lastband_number = 0;
        int count_bands = 0;
        for (int i = 0; i < dsize; i++)
        {
            if (i < nbands % dsize)
            {
                // mohan modify 2010-07-05
                pv->loc_sizes[i] = nbands / dsize + 1;
            }
            else
            {
                pv->loc_sizes[i] = nbands / dsize;
            }
            count_bands += pv->loc_sizes[i];
            if (count_bands >= nbands)
            {
                pv->lastband_in_proc = i;
                pv->lastband_number = nbands - (count_bands - pv->loc_sizes[i]);
                break;
            }
        }
    }
    else
    {
        pv->loc_size = nlocal / dsize;

        // mohan add 2012-03-29
        if (pv->loc_size == 0)
        {
            ofs_warning << " loc_size=0" << " in proc " << myrank + 1 << std::endl;
            ModuleBase::WARNING_QUIT("ORB_control::set_parameters", "nbands < ncpus");
        }

        if (drank < nlocal % dsize)
        {
            pv->loc_size += 1;
        }
        if (pv->testpb) ModuleBase::GlobalFunc::OUT(ofs_running, "local size", pv->loc_size);

        // set loc_sizes
        delete[] pv->loc_sizes;
        pv->loc_sizes = new int[dsize];
        ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, dsize);

        pv->lastband_in_proc = 0;
        pv->lastband_number = 0;
        int count_bands = 0;
        for (int i = 0; i < dsize; i++)
        {
            if (i < nlocal % dsize)
            {
                // mohan modify 2010-07-05
                pv->loc_sizes[i] = nlocal / dsize + 1;
            }
            else
            {
                pv->loc_sizes[i] = nlocal / dsize;
            }
            count_bands += pv->loc_sizes[i];
            if (count_bands >= nbands)
            {
                pv->lastband_in_proc = i;
                pv->lastband_number = nbands - (count_bands - pv->loc_sizes[i]);
                break;
            }
        }
    }//xiaohui add 2014-12-21

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "lastband_in_proc", pv->lastband_in_proc);
    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "lastband_number", pv->lastband_number);

    return;
}


#ifdef __MPI
// creat the 'comm_2D' stratege.
void ORB_control::mpi_creat_cart(MPI_Comm* comm_2D,
    int prow, int pcol, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ORB_control", "mpi_creat_cart");
    // the matrix is divided as ( dim[0] * dim[1] )
    int dim[2];
    int period[2] = { 1,1 };
    int reorder = 0;
    dim[0] = prow;
    dim[1] = pcol;

    if (this->ParaV.testpb) ofs_running << " dim = " << dim[0] << " * " << dim[1] << std::endl;

    MPI_Cart_create(DIAG_WORLD, 2, dim, period, reorder, comm_2D);
    return;
}
#endif

#ifdef __MPI
int ORB_control::mat_2d(MPI_Comm vu,
    const int& M_A,
    const int& N_A,
    const int& nb,
    LocalMatrix& LM,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "mat_2d");

    Parallel_Orbitals* pv = &this->ParaV;

    int dim[2];
    int period[2];
    int coord[2];
    int i, j, k, end_id;
    int block;

    // (0) every processor get it's id on the 2D comm
    // : ( coord[0], coord[1] )
    MPI_Cart_get(vu, 2, dim, period, coord);

    // (1.1) how many blocks at least
    // eg. M_A = 6400, nb = 64;
    // so block = 10;
    block = M_A / nb;

    // (1.2) If data remain, add 1.
    if (block * nb < M_A)
    {
        block++;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Total Row Blocks Number", block);

    // mohan add 2010-09-12
    if (dim[0] > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of row blocks is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ORB_control::mat_2d", "some processor has no row blocks, try a smaller 'nb2d' parameter.");
        }
    }

    // (2.1) row_b : how many blocks for this processor. (at least)
    LM.row_b = block / dim[0];

    // (2.2) row_b : how many blocks in this processor.
    // if there are blocks remain, some processors add 1.
    if (coord[0] < block % dim[0])
    {
        LM.row_b++;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local Row Block Number", LM.row_b);

    // (3) end_id indicates the last block belong to
    // which processor.
    if (block % dim[0] == 0)
    {
        end_id = dim[0] - 1;
    }
    else
    {
        end_id = block % dim[0] - 1;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Ending Row Block in processor", end_id);

    // (4) row_num : how many rows in this processors :
    // the one owns the last block is different.
    if (coord[0] == end_id)
    {
        LM.row_num = (LM.row_b - 1) * nb + (M_A - (block - 1) * nb);
    }
    else
    {
        LM.row_num = LM.row_b * nb;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local rows (including nb)", LM.row_num);

    // (5) row_set, it's a global index :
    // save explicitly : every row in this processor
    // belongs to which row in the global matrix.
    LM.row_set.resize(LM.row_num);
    j = 0;
    for (i = 0; i < LM.row_b; i++)
    {
        for (k = 0; k < nb && (coord[0] * nb + i * nb * dim[0] + k < M_A); k++, j++)
        {
            LM.row_set[j] = coord[0] * nb + i * nb * dim[0] + k;
            // ofs_running << " j=" << j << " row_set=" << LM.row_set[j] << std::endl;
        }
    }

    // the same procedures for columns.
    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Total Col Blocks Number", block);

    if (dim[1] > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of column blocks is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ORB_control::mat_2d", "some processor has no column blocks.");
        }
    }

    LM.col_b = block / dim[1];
    if (coord[1] < block % dim[1])
    {
        LM.col_b++;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local Row Block Number", LM.col_b);

    if (block % dim[1] == 0)
    {
        end_id = dim[1] - 1;
    }
    else
    {
        end_id = block % dim[1] - 1;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Ending Row Block in processor", end_id);

    if (coord[1] == end_id)
    {
        LM.col_num = (LM.col_b - 1) * nb + (M_A - (block - 1) * nb);
    }
    else
    {
        LM.col_num = LM.col_b * nb;
    }

    if (pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running, "Local columns (including nb)", LM.row_num);

    LM.col_set.resize(LM.col_num);

    j = 0;
    for (i = 0; i < LM.col_b; i++)
    {
        for (k = 0; k < nb && (coord[1] * nb + i * nb * dim[1] + k < M_A); k++, j++)
        {
            LM.col_set[j] = coord[1] * nb + i * nb * dim[1] + k;
        }
    }
    LM.col_pos = 0;
    LM.row_pos = 0;

    // for wavefuncton , calculate nbands_loc
    block = N_A / nb;
    if (block * nb < N_A)
    {
        block++;
    }
    if (dim[1] > block)
    {
        ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
        ofs_warning << " but, the number of bands-row-block is " << block << std::endl;
        if (nb > 1)
        {
            return 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ORB_control::mat_2d", "some processor has no bands-row-blocks.");
        }
    }
    int col_b_bands = block / dim[1];
    if (coord[1] < block % dim[1])
    {
        col_b_bands++;
    }
    if (block % dim[1] == 0)
    {
        end_id = dim[1] - 1;
    }
    else
    {
        end_id = block % dim[1] - 1;
    }
    if (coord[1] == end_id)
    {
        pv->ncol_bands = (col_b_bands - 1) * nb + (N_A - (block - 1) * nb);
    }
    else
    {
        pv->ncol_bands = col_b_bands * nb;
    }
    pv->nloc_wfc = pv->ncol_bands * LM.row_num;

    pv->nloc_Eij= pv->ncol_bands * pv->ncol_bands;
    
    return 0;
}
#endif
