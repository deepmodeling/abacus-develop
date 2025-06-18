#include <gtest/gtest.h>
#include <complex>
#include <memory>
#include <vector>
#include <cmath>

#include "module_io/cal_pLpR.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_basis/module_nao/radial_collection.h"
#include "source_base/spherical_bessel_transformer.h"
#include "source_base/ylm.h"

#define DOUBLETHRESHOLD 1e-12

class CalpLpRTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        ModuleBase::SphericalBesselTransformer sbt(true);

        this->orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        this->orb_->build(1, &forb_, 'o');
        this->orb_->set_transformer(sbt);

        const double rcut = 8.0;
        const int ngrid = int(rcut / 0.01) + 1;
        const double cutoff = 2.0 * rcut;
        this->orb_->set_uniform_grid(true, ngrid, cutoff, 'i', true);

        this->calculator_ = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        this->calculator_->tabulate(*orb_, *orb_, 'S', ngrid, cutoff);

        // Initialize Ylm coefficients
        ModuleBase::Ylm::set_coefficients();
    }

    void TearDown() override
    {
        // Cleanup code here, if needed
    }

    const std::string forb_ = "../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";

    std::unique_ptr<TwoCenterIntegrator> calculator_;
    std::unique_ptr<RadialCollection> orb_;
};

TEST_F(CalpLpRTest, CalLzijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, mi = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, mj = 0; // self, the first s
    
    // <s|Lz|s> = 0: no magnetic moment
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <s|Lz|p> = 0: orthogonal
    jl = 1;
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <p(m=-1)|Lz|p(m=0)> = 0: orthogonal
    il = 1; mi = -1; jl = 1; mj = 0;
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <p(m=1)|Lz|p(m=1)> = 1: same
    il = 1; mi = 1; jl = 1; mj = 1;
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <d(m=-1)|Lz|d(m=0)> = 0: orthogonal
    il = 2; mi = -1; jl = 2; mj = 0;
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <d(m=1)|Lz|d(m=1)> = 1: same
    il = 2; mi = 1; jl = 2; mj = 1;
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, mi, jt, ja, jl, jz, mj, vR);
    EXPECT_NEAR(out.real(), 1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
}

TEST_F(CalpLpRTest, CalLxijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0; // self, the first s

    // <s|Lx|s> = 0: no anisotropy
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <s|Lx|p> = 0: orthogonal
    jl = 1;
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <p(m=-1)|Lx|p(m=0)> = 0.5: Lx|p(m=0)> = 1/2 (|p(m=-1)> + |p(m=1)>)
    il = 1; im = -1; jl = 1; jm = 0;
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.5*sqrt(2.0), DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <p(m=1)|Lx|p(m=1)> = 0: expectation value is 0
    il = 1; im = 1; jl = 1; jm = 1;
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <d(m=-1)|Lx|d(m=0)> = 0.5: Lx|d(m=0)> = 1/2 (|d(m=-1)> + |d(m=1)>)
    il = 2; im = -1; jl = 2; jm = 0;
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.5*sqrt(6.0), DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
}

TEST_F(CalpLpRTest, CalLyijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0; // self, the first s

    // <s|Ly|s> = 0: no anisotropy
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <s|Ly|p> = 0: orthogonal
    jl = 1;
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <p(m=-1)|Ly|p(m=0)> = -i/2: Ly|p(m=0)> = -i/2 (|p(m=1)> - |p(m=-1)>)
    il = 1; im = -1; jl = 1; jm = 0;
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.5*sqrt(2.0), DOUBLETHRESHOLD);
    // <p(m=1)|Ly|p(m=1)> = 0: expectation value is 0
    il = 1; im = 1; jl = 1; jm = 1;
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    // <d(m=-1)|Ly|d(m=0)> = -i/2: Ly|d(m=0)> = -i/2 (|d(m=1)> - |d(m=-1)>)
    il = 2; im = -1; jl = 2; jm = 0;
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.5*sqrt(6.0), DOUBLETHRESHOLD);
}

int main(int argc, char **argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    return 0;
}