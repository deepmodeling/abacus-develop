/**
 * Unit-test of cal_pLpR.cpp
 * 
 * The representation matrices for Lx, Ly and Lz operators (SO(3) generators).
 * On basis {s}:
 * Lx: Ly: Lz:
 * 0   0   0
 * 
 * On basis {px, py, pz}:
 * Lx:             Ly:             Lz:
 *    px  py  pz      px  py  pz      px  py  pz
 * px  0           px  0       i   px  0  -i
 * py      0  -i   py      0       py  i   0
 * pz      i   0   pz -i       0   pz          0
 * 
 * will test if the calculated values are close to the above results.
 */

#include <gtest/gtest.h>
#include <complex>
#include <memory>
#include <vector>
#include <cmath>

#include "source_io/module_hs/cal_pLpR.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_basis/module_nao/radial_collection.h"
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
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0;
    
    // l=0
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);

    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[1] = {0.0, -1.0}; ans[3] = {0.0, 1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (auto im_: m) {
        for (auto jm_: m) {
            im = im_; jm = jm_;
            out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), DOUBLETHRESHOLD);
            idx++;
        }
    }
}

TEST_F(CalpLpRTest, CalLxijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0;

    // l=0
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);

    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[5] = {0.0, 1.0}; ans[7] = {0.0, -1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (auto im_: m) {
        for (auto jm_: m) {
            im = im_; jm = jm_;
            out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), DOUBLETHRESHOLD);
            idx++;
        }
    }
}

TEST_F(CalpLpRTest, CalLyijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0; // self, the first s

    // l=0
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out.imag(), 0.0, DOUBLETHRESHOLD);
    
    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[2] = {0.0, -1.0}; ans[6] = {0.0, 1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (auto im_: m) {
        for (auto jm_: m) {
            im = im_; jm = jm_;
            out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), DOUBLETHRESHOLD);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), DOUBLETHRESHOLD);
            idx++;
        }
    }
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