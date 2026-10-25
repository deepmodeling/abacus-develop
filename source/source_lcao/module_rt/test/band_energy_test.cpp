#include "source_lcao/module_rt/band_energy.h"

#include "tddft_test.h"

#include <cstdio>
#include <fstream>
#include <iterator>

namespace
{
class BandEnergyTest : public testing::Test
{
  protected:
    using Complex = std::complex<double>;
    const int nlocal = 4;
    const int nband = 3;
    Parallel_Orbitals pv;
    ct::Tensor h{ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({16})};
    ct::Tensor psi{ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({12})};
    ct::Tensor energies{ct::DataType::DT_DOUBLE, ct::TensorShape({3})};
    const module_rt::EijOutputInfo output_info{0.5, 2, 3, 1, 1, 1};
    std::string output_path;
    std::ofstream output;

    void SetUp() override
    {
        // Reuse the existing test executable's 1x1 BLACS grid.
        ASSERT_EQ(nprow, 1);
        ASSERT_EQ(npcol, 1);
        ASSERT_EQ(pv.set(nlocal, nlocal, 1, ictxt), 0);
        pv.ncol_bands = nband;
        pv.nrow_bands = nband;
        pv.nloc_wfc = nlocal * nband;
        pv.set_desc_wfc_Eij(nlocal, nband, nlocal);

        // Identity columns make Eij the leading 3x3 block of this column-major H.
        // The fourth basis function also exercises the LAPACK output stride.
        h.zero();
        psi.zero();
        const double diagonal[]{3.0, 8.0, 10.0, 12.0};
        Complex* matrix = h.data<Complex>();
        for (int i = 0; i < nlocal; ++i)
        {
            matrix[i * nlocal + i] = diagonal[i];
        }
        matrix[4] = Complex(0.5, 0.125);
        matrix[8] = Complex(0.125, 0.75);
        matrix[9] = Complex(0.125, 0.25);
        matrix[1] = std::conj(matrix[4]);
        matrix[2] = std::conj(matrix[8]);
        matrix[6] = std::conj(matrix[9]);
        for (int i = 0; i < nband; ++i)
        {
            psi.data<Complex>()[i * nlocal + i] = 1.0;
        }

        const std::string name = testing::UnitTest::GetInstance()->current_test_info()->name();
        output_path = testing::TempDir() + "band_energy_" + name + ".txt";
        output.open(output_path);
        ASSERT_TRUE(output.is_open());
    }

    void TearDown() override
    {
        output.close();
        std::remove(output_path.c_str());
    }

    void check_result()
    {
        EXPECT_NEAR(energies.data<double>()[0], 3.0, 1e-12);
        EXPECT_NEAR(energies.data<double>()[1], 8.0, 1e-12);
        EXPECT_NEAR(energies.data<double>()[2], 10.0, 1e-12);
        output.close();
        std::ifstream input(output_path);
        const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
        // Both conjugate entries must keep the sub-threshold component; equality passes the threshold.
        EXPECT_NE(text.find("i =  1, j =  2, Eij =   5.00000000e-01 + 1.25000000e-01 i"), std::string::npos);
        EXPECT_NE(text.find("i =  2, j =  1, Eij =   5.00000000e-01 - 1.25000000e-01 i"), std::string::npos);
        EXPECT_NE(text.find("i =  1, j =  3, Eij =   1.25000000e-01 + 7.50000000e-01 i"), std::string::npos);
        EXPECT_NE(text.find("i =  3, j =  1, Eij =   1.25000000e-01 - 7.50000000e-01 i"), std::string::npos);
        EXPECT_EQ(text.find("i =  2, j =  3"), std::string::npos);
        EXPECT_EQ(text.find("i =  3, j =  2"), std::string::npos);
    }
};

TEST_F(BandEnergyTest, ComplexEijScalapack)
{
    module_rt::compute_ekb(&pv, nband, nlocal, h.data<Complex>(), psi.data<Complex>(), energies.data<double>(), output, output_info);
    check_result();
}

TEST_F(BandEnergyTest, ComplexEijLapack)
{
    module_rt::compute_ekb_tensor_lapack<base_device::DEVICE_CPU>(&pv, nband, nlocal, h, psi, energies, output, output_info);
    check_result();
}
} // namespace
