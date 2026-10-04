#ifndef SCCS_TEST_H
#define SCCS_TEST_H

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/pw_basis.h"
#include <gtest/gtest.h>

#include <cmath>

namespace SccsTest
{
extern int pool_size;
extern int pool_rank;

class PwTest : public testing::Test
{
protected:
    PwTest() : basis("cpu", "double") {}

    void SetUp() override
    {
#ifdef __MPI
        basis.initmpi(pool_size, pool_rank, POOL_WORLD);
#endif
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0);
        basis.initgrids(length, lattice, 20.0);
        basis.initparameters(false, 20.0, 1, false);
        basis.setuptransform();
        basis.collect_local_pw();
        tpiba = ModuleBase::TWO_PI / length;
    }

    // A cosine along x or y, in the PW local real-grid ordering.
    std::vector<double> cosine_mode(int direction) const
    {
        std::vector<double> values(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const int ix = ir / (basis.ny * basis.nplane);
            const int iy = ir / basis.nplane - ix * basis.ny;
            const int coordinate = direction == 0 ? ix : iy;
            const int count = direction == 0 ? basis.nx : basis.ny;
            const double angle = ModuleBase::TWO_PI * coordinate / count;
            values[ir] = std::cos(angle);
        }
        return values;
    }

    const double length = 10.0;
    double tpiba = 0.0;
    ModulePW::PW_Basis basis;
    std::string error;
};
} // namespace SccsTest

#endif
