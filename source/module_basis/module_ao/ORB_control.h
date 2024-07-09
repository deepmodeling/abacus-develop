#ifndef ORB_CONTROL_H
#define ORB_CONTROL_H

#include "ORB_read.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"
#include "parallel_orbitals.h"

class ORB_control
{
  public:
    /// use this when need to init 2D-division-info (ParaV)
    /// which is used next in wfc/matrix
    ORB_control(const bool& gamma_only_in,
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

    /// set 2D-block-cyclic division according to the basis
    void setup_2d_division(std::ofstream& ofs_running, std::ofstream& ofs_warning);

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


};
#endif
