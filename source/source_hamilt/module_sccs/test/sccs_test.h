#ifndef SCCS_TEST_H
#define SCCS_TEST_H

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_comm.h"
#include "source_basis/module_pw/pw_basis.h"
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace SccsTest
{
// Model charges of the open-boundary tests, mean + amplitude * cosine mode.
// The mean gives a monopole and the mode has an analytic potential; the
// operators are linear, so only the shape matters, and the values are of
// valence-density size with the density kept positive.
const double charge_mean = 1e-3;
const double charge_amplitude = 2e-4;
// A second charge, not proportional to the first, for symmetry checks.
const double probe_mean = 3e-4;
const double probe_amplitude = 1e-4;
// Cavity density inside the water-neutral switching range (1e-4 to 5e-3), so
// the dielectric varies in space; the small amplitude keeps it mildly
// nonuniform and the sqrt-CG short.
const double cavity_mean = 1e-3;
const double cavity_amplitude = 2e-5;

extern int pool_size;
extern int pool_rank;
// Probe kernel of a cavity without the solvent-aware filling.
const std::vector<double> no_probe;

class PwTest : public testing::Test
{
protected:
    PwTest() : basis("cpu", "double") {}

    void SetUp() override
    {
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0);
        set_up_basis(lattice);
    }

    // Rows of lattice in units of length; ecut 20 Ry.
    void set_up_basis(const ModuleBase::Matrix3& lattice)
    {
#ifdef __MPI
        basis.initmpi(pool_size, pool_rank, POOL_WORLD);
#endif
        basis.initgrids(length, lattice, 20.0);
        finish_basis();
    }

    // Explicit FFT grid. A pool with more ranks than nz leaves ranks without
    // real-space grid points (nrxx == 0); they still own G vectors.
    void set_up_grid(const ModuleBase::Matrix3& lattice, int nx, int ny, int nz)
    {
#ifdef __MPI
        basis.initmpi(pool_size, pool_rank, POOL_WORLD);
#endif
        basis.initgrids(length, lattice, nx, ny, nz);
        finish_basis();
    }

    void finish_basis()
    {
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
// A cell short along c on a 12 x 12 x 2 grid: with three pool ranks the
// third owns no z plane, while the G sphere keeps sticks for every rank.
class EmptyPartitionTest : public PwTest
{
protected:
    void SetUp() override
    {
        set_up_grid(lattice, 12, 12, 2);
    }

    const ModuleBase::Matrix3 lattice = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                                            0.0, 1.0, 0.0,
                                                            0.0, 0.0, 0.25);
};
} // namespace SccsTest

#endif
