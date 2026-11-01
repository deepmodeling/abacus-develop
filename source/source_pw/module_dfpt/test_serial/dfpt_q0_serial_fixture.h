#ifndef DFPT_Q0_SERIAL_FIXTURE_H
#define DFPT_Q0_SERIAL_FIXTURE_H

// Shared serial-side fixture for the DFPT_Q0 test files split along the
// source boundaries (test_dfpt_q0.cpp, test_dfpt_q0_pos.cpp). Runs
// without __MPI on the DFPTSerialBase shared FFT grid like the other
// DFPT serial tests; all references are closed-form or operator finite
// differences, no ground-state solver is involved.

#include "dfpt_serial_fixture.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include "source_pw/module_dfpt/dfpt_pert.h"
#include "source_pw/module_dfpt/dfpt_q0.h"
#include "source_cell/stru_fac.h"

class DFPTQ0SerialTest : public DFPTSerialBase
{
  protected:
    Structure_Factor sf_;
    ModuleDFPT::DFPT_Pert pert_;
    ModuleDFPT::DFPT_Q0 q0_;

    // this fixture is Gamma-only: shadow the generic default k of the base
    const ModuleBase::Vector3<double> k_d_{0.0, 0.0, 0.0};
    const ModuleBase::Vector3<double> gx_{0.1, 0.0, 0.0}; // 1/lat0 units
    const ModuleBase::Vector3<double> gy_{0.0, 0.1, 0.0};

    void SetUp() override
    {
        SetUpCell();
        SetupBases(k_d_, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), 4);
        pert_.init(ucell_, &pw_rho_, &pw_wfc_, sf_);
        q0_.init(ucell_, &pw_rho_, &pw_wfc_, &pert_);
    }

    // wfc-basis index of the reciprocal vector (ix, iy, iz)/a at Gamma
    int IgOf(int ix, int iy, int iz) const
    {
        const int npwk = pw_wfc_.npwk[0];
        for (int ig = 0; ig < npwk; ++ig)
        {
            const ModuleBase::Vector3<double> g = pw_wfc_.getgpluskcar(0, ig);
            if (std::llround(g.x * a_) == ix && std::llround(g.y * a_) == iy && std::llround(g.z * a_) == iz)
            {
                return ig;
            }
        }
        return -1;
    }
};

#endif // DFPT_Q0_SERIAL_FIXTURE_H
