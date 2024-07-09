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
    ORB_control(const int& nlocal_in,
                const int& nbands_in,
                const int& nb2d_in);

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

    // -------------------------------------------------------------------------
    // note: ORB_control orb_con is now a member in ESolver_KS_LCAO
    // ("friend class ESolver_KS_LCAO;" will cause not-defined problem).
    // These variables is set in in ESolver_KS_LCAO
    // and can only be visited in ESolver_KS_LCAO.
    // -------------------------------------------------------------------------
    int nlocal = 0;
    int nbands = 0;
    int nb2d = 0;
};
#endif
