#ifndef ORB_CONTROL_H 
#define ORB_CONTROL_H 

#include "module_io/input.h"
#include "module_cell/unitcell.h"
#include "parallel_orbitals.h"
#include "ORB_gen_tables.h"
#include "ORB_read.h"

class ORB_control
{
public:

    /// use this when need to init 2D-division-info (ParaV)
    /// which is used next in wfc/matrix
    ORB_control(
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
        const std::string& ks_solver_in);

    /// use this when only need to calculate 
    /// 2-center-integral of read orbitals
    ORB_control();

    ~ORB_control();

    void Init(Input& inp, UnitCell& ucell);

    //first step: read orbital file
    void read_orb_first(
        std::ofstream& ofs_in,
        LCAO_Orbitals& orb,
        const int& ntype, // mohan add 2021-04-26
        const std::string& orbital_dir,  // liuyu add 2023-04-06
        const std::string* orbital_file,  // liuyu add 2023-04-06
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
    );
    /// Generate the S(overlap),T,NL matrix.
    void set_orb_tables(
        std::ofstream& ofs_in,
        ORB_gen_tables& OGT,
        LCAO_Orbitals& orb,
        const double& lat0,
        const bool& deepks_setorb,
        const int& Lmax_exx,
        const int& nprojmax,
        const int* nproj,
        const Numerical_Nonlocal* beta_);

    void clear_after_ions(
        ORB_gen_tables& OGT,
        LCAO_Orbitals& orb,
        const bool& deepks_setorb,
        const int* nproj_);

    ///set 2D-block-cyclic division according to the basis
    void setup_2d_division(std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

    // #ifdef __MPI
    //     void readin(const std::string& fa, const std::string& fb, const int& nlocal, double* eigen, double* eigvr);
    // #endif

    Parallel_Orbitals ParaV;

    bool setup_2d = false;

    // -------------------------------------------------------------------------
    // note: ORB_control orb_con is now a member in ESolver_KS_LCAO
    // ("friend class ESolver_KS_LCAO;" will cause not-defined problem).
    // These variables is set in in ESolver_KS_LCAO 
    // and can only be visited in ESolver_KS_LCAO.
    // -------------------------------------------------------------------------
    bool gamma_only = 1;
    int nlocal = 0;
    int nbands = 0;
    int nspin = 1;
    int dsize = 0;
    int nb2d = 0;
    int dcolor = 0;
    int drank = 0;
    int myrank = 0;
    std::string calculation;
    std::string ks_solver;
    std::ofstream ofs_running;
    std::ofstream ofs_warning;

private:
    void divide_HS_2d
    (
#ifdef __MPI
        MPI_Comm DIAG_WORLD,
#endif
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

#ifdef __MPI
    int mpi_comm_rows, mpi_comm_cols;
#endif

    void set_parameters(std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

    void set_trace(std::ofstream& ofs_running);

#ifdef __MPI
    int mat_2d(MPI_Comm vu,
        const int& M_A, const int& N_A,
        const int& NB, LocalMatrix& loc_A,
        std::ofstream& ofs_running,
        std::ofstream& ofs_warning);

#endif

};
#endif
