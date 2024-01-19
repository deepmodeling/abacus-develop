#include "ORB_control.h"
#include "ORB_gen_tables.h"
#include "module_base/timer.h"
#include "module_base/parallel_common.h"
#include "module_base/lapack_connector.h"
#include "module_base/blacs_connector.h"
#include "module_base/memory.h"
#include "module_base/parallel_global.h"

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
{
    Cblacs_exit(1); //delete global variables in cblacs but do not close MPI
}

void ORB_control::read_orb_first(
    std::ofstream& ofs_in,
    LCAO_Orbitals& orb,
    const int& ntype, // mohan add 2021-04-26
    const std::string& orbital_dir,  // liuyu add 2023-04-06
    const std::string *orbital_file,  // liuyu add 2023-04-06
    const std::string& descriptor_file,  // liuyu add 2023-04-06
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

    // liuyu add 2023-04-06
    if (my_rank == 0 && !orb.read_in_flag)
    {
        orb.read_in_flag = true;
        orb.descriptor_file = descriptor_file;
        for(int it = 0; it < ntype; ++it)
        {
            std::string ofile = orbital_dir + orbital_file[it];
            orb.orbital_file.push_back(ofile);
        }
    }
#ifdef __MPI
    orb.bcast_files(ntype, my_rank);
#endif

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
    
    /// generate overlap & kinetic table
    OGT.gen_tables(ofs_in, orb, Lmax_exx, deepks_setorb, nprojmax, nproj, beta_);
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
    if (ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver" || ks_solver == "cg_in_lcao")
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

    // determine whether 2d-division or not according to ks_solver
    bool div_2d;
    if (ks_solver == "lapack" || ks_solver == "cg" || ks_solver == "dav") div_2d = false;
#ifdef __MPI
    else if (ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver" || ks_solver == "cg_in_lcao") div_2d = true;
#endif
    else
    {
        std::cout << " Parallel Orbial, DIAGO_TYPE = " << ks_solver << std::endl;
        ModuleBase::WARNING_QUIT("Parallel_Orbitals::set_global2local", "Check ks_solver.");
    }
    // (2) set the trace, then we can calculate the nnr.
    // for 2d: calculate po.nloc first, then global2local_row and global2local_col
    // for O(N): calculate the three together.
    this->ParaV.set_global2local(nlocal, nlocal, div_2d, ofs_running);
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

// divide the H&S matrix using 2D block algorithms.
void ORB_control::divide_HS_2d(
#ifdef __MPI
    MPI_Comm DIAG_WORLD,
#endif
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "divide_HS_2d");
    assert(nlocal > 0);
    assert(dsize > 0);
    Parallel_Orbitals* pv = &this->ParaV;

    if (dcolor != 0)
        return; // mohan add 2012-01-13

    // get the 2D index of computer.
    pv->dim0 = (int)sqrt((double)dsize); // mohan update 2012/01/13
    // while (GlobalV::NPROC_IN_POOL%dim0!=0)

    if (ks_solver == "cusolver")
    {
        pv->dim0 = 1; pv->dim1 = dsize;
    } // Xu Shu add 2022-03-25
    else
        pv->set_proc_dim(dsize);

    if (pv->testpb)
        ModuleBase::GlobalFunc::OUT(ofs_running, "dim0", pv->dim0);
    if (pv->testpb)
        ModuleBase::GlobalFunc::OUT(ofs_running, "dim1", pv->dim1);

#ifdef __MPI
    // mohan add 2011-04-16
#ifdef __DEBUG
assert(nb2d > 0);
#endif
    pv->set_block_size(nb2d); // mohan add 2010-06-28

    if (ks_solver == "cusolver")
        pv->set_block_size(1); // Xu Shu add 2022-03-25
    ModuleBase::GlobalFunc::OUT(ofs_running, "nb2d", pv->get_block_size());

    this->set_parameters(ofs_running, ofs_warning);

    // call mpi_creat_cart
    pv->mpi_create_cart(DIAG_WORLD);

    int try_nb = pv->set_local2global(nlocal, nlocal, ofs_running, ofs_warning);
    try_nb = pv->set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);
    if (try_nb == 1)
    {
        ofs_running << " parameter nb2d is too large: nb2d = " << pv->get_block_size() << std::endl;
        ofs_running << " reset nb2d to value 1, this set would make the program keep working but maybe get slower "
                       "during diagonalization."
                    << std::endl;
        pv->set_block_size(1);
        try_nb = pv->set_local2global(nlocal, nlocal, ofs_running, ofs_warning);
        try_nb = pv->set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);
    }

    // init blacs context for genelpa
    if (ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver" || ks_solver == "cg_in_lcao")
    {
        pv->set_desc(nlocal, nlocal, pv->nrow);
        pv->set_desc_wfc_Eij(nlocal, nbands, pv->nrow);
    }
#else // single processor used.
    pv->set_serial(nlocal, nlocal);
    this->set_parameters(ofs_running, ofs_warning);
#endif

    assert(pv->nloc > 0);
    if (pv->testpb)
        ModuleBase::GlobalFunc::OUT(ofs_running, "this->nrow", pv->nrow);
    if (pv->testpb)
        ModuleBase::GlobalFunc::OUT(ofs_running, "this->ncol", pv->ncol);
    if (pv->testpb)
        ModuleBase::GlobalFunc::OUT(ofs_running, "nloc", pv->nloc);
    return;
}
