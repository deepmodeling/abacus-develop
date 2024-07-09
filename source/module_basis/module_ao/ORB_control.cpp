#include "ORB_control.h"

#include "module_base/blacs_connector.h"
#include "module_base/lapack_connector.h"
#include "module_base/memory.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_global.h"
#include "module_base/timer.h"

//#include "build_st_pw.h"

ORB_control::ORB_control(const bool& gamma_only_in,
                         const int& nlocal_in,
                         const int& nbands_in,
                         const int& nspin_in,
                         const int& dsize_in,
                         const int& nb2d_in,
                         const int& dcolor_in,
                         const int& drank_in,
                         const int& myrank_in,
                         const std::string& calculation_in,
                         const std::string& ks_solver_in)
    : gamma_only(gamma_only_in), nlocal(nlocal_in), nbands(nbands_in), nspin(nspin_in), dsize(dsize_in), nb2d(nb2d_in),
      dcolor(dcolor_in), drank(drank_in), myrank(myrank_in), calculation(calculation_in), ks_solver(ks_solver_in),
      setup_2d(true)
{
    this->ParaV.nspin = nspin_in;
}

ORB_control::ORB_control() : setup_2d(false)
{
}
ORB_control::~ORB_control()
{
#ifdef __MPI
    Cblacs_exit(1); // delete global variables in cblacs but do not close MPI
#endif
}


void ORB_control::setup_2d_division(std::ofstream& ofs_running, std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "setup_2d_division");
    ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << std::endl;

#ifdef __MPI
    // storage form of H and S matrices on each processor
    // is determined in 'divide_HS_2d' subroutine

    ModuleBase::GlobalFunc::OUT(ofs_running, "nb2d", ParaV.get_block_size());
    int try_nb = ParaV.init(nlocal, nlocal, nb2d, DIAG_WORLD);
    try_nb += ParaV.set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);
    if (try_nb != 0)
    {
        ofs_running << " parameter nb2d is too large: nb2d = " << ParaV.get_block_size() << std::endl;
        ofs_running << " reset nb2d to value 1, this set would make the program keep working but maybe get slower "
                       "during diagonalization."
                    << std::endl;

        ParaV.set(nlocal, nlocal, 1, ParaV.comm_2D, ParaV.blacs_ctxt);
        try_nb = ParaV.set_nloc_wfc_Eij(nbands, ofs_running, ofs_warning);
    }

    // init blacs context for genelpa
    if (ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver" || ks_solver == "cusolvermp"
        || ks_solver == "cg_in_lcao" || ks_solver == "pexsi")
    {
        ParaV.set_desc_wfc_Eij(nlocal, nbands, ParaV.nrow);
    }

#else
    ParaV->set_serial(nlocal, nlocal);
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06
#endif
}

