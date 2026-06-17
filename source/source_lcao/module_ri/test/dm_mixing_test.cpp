#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_lcao/module_ri/Mix_DMk_2D.h"

/************************************************
 *  unit test of charge_mixing.cpp & Mix_DMk_2D.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Mix_DMk_2D::mix:
 *      mix the density matrix data according to the set mixing mode.
 *
 */

class DM_Mixing_Test : public ::testing::Test
{
  public:
    DM_Mixing_Test()
    {
        mixing = new Base_Mixing::Broyden_Mixing(ndim, mixing_beta);
        mix_data_vector = std::vector<std::vector<double>>(2);
        mix_complexdata_vector = std::vector<std::vector<std::complex<double>>>(3);
        mix_data_vector[0].resize(nr * nc);
        mix_data_vector[1].resize(nr * nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_data_vector[0][i * nc + j] = i * nc + j;
                mix_data_vector[1][i * nc + j] = i * nc + j + 0.2;
            }
        }
        mix_complexdata_vector[0].resize(nr * nc);
        mix_complexdata_vector[1].resize(nr * nc);
        mix_complexdata_vector[2].resize(nr * nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_complexdata_vector[0][i * nc + j] = std::complex<double>{ double(i), double(j) };
                mix_complexdata_vector[1][i * nc + j] = std::complex<double>{ double(i), double(j) + 0.2 };
                mix_complexdata_vector[2][i * nc + j] = std::complex<double>{ double(i) + 0.8, double(j) };
            }
        }
    };
    ~DM_Mixing_Test()
    {
        delete mixing;
    };
    Base_Mixing::Mixing* mixing = nullptr;
    const int nr = 2;
    const int nc = 2;
    const int ndim = 1;
    const double mixing_beta = 0.3;

  protected:
    std::vector<std::vector<double>> mix_data_vector;
    std::vector<std::vector<std::complex<double>>> mix_complexdata_vector;
};

TEST_F(DM_Mixing_Test, Mix_DMk_2D)
{
    //Gamma only
    Mix_DMk_2D<double> mix_dmk_gamma;
    mix_dmk_gamma.set_nks(1);
    mix_dmk_gamma.set_mixing_plain(1.0);
    std::vector<std::vector<std::vector<double>>> dm_gamma(2);
    dm_gamma[0] = std::vector<std::vector<double>>(1);
    dm_gamma[0][0] = mix_data_vector[0];
    dm_gamma[1] = std::vector<std::vector<double>>(1);
    dm_gamma[1][0] = mix_data_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk_gamma.mix(dm_gamma[istep], (istep == 0));
    }
    std::vector<const std::vector<double>*> dm_gamma_out = mix_dmk_gamma.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(dm_gamma_out[0][0][i * nc + j], mix_data_vector[1][i * nc + j]);
        }
    }

    // not Gamma only
    Mix_DMk_2D<std::complex<double>> mix_dmk;
    mix_dmk.set_nks(1);
    mix_dmk.set_mixing_plain(1.0);
    std::vector<std::vector<std::vector<std::complex<double>>>> dm(2);
    dm[0] = std::vector<std::vector<std::complex<double>>>(1);
    dm[0][0] = mix_complexdata_vector[0];
    dm[1] = std::vector<std::vector<std::complex<double>>>(1);
    dm[1][0] = mix_complexdata_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk.mix(dm[istep], (istep == 0));
    }
    std::vector<const std::vector<std::complex<double>>*> dm_out = mix_dmk.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].real(), mix_complexdata_vector[1][i * nc + j].real());
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].imag(), mix_complexdata_vector[1][i * nc + j].imag());
        }
    }

    // Shared Broyden mix
    Mix_DMk_2D<std::complex<double>> mix_dmk_broyden;
    mix_dmk_broyden.set_nks(1);
    mix_dmk_broyden.set_mixing(mixing);
    mixing->coef = { 1.1, -0.1 };
    std::vector<std::vector<std::vector<std::complex<double>>>> dm_broyden(3);
    for (int istep = 0; istep < 3; ++istep)
    {
        dm_broyden[istep] = std::vector<std::vector<std::complex<double>>>(1);
        dm_broyden[istep][0] = mix_complexdata_vector[istep];
        mix_dmk_broyden.mix(dm_broyden[istep], (istep == 0));
    }
    std::vector<const std::vector<std::complex<double>>*> dm_broyden_out = mix_dmk_broyden.get_DMk_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            std::complex<double> first_step_result
                = (1 - mixing_beta) * mix_complexdata_vector[0][i * nc + j]
                  + mixing_beta * mix_complexdata_vector[1][i * nc + j];
            std::complex<double> second_step_result
                = (1 - mixing_beta) * first_step_result + mixing_beta * mix_complexdata_vector[2][i * nc + j];
            std::complex<double> ref = second_step_result * mixing->coef[1] + first_step_result * mixing->coef[0];
            EXPECT_DOUBLE_EQ(dm_broyden_out[0][0][i * nc + j].real(), ref.real());
            EXPECT_DOUBLE_EQ(dm_broyden_out[0][0][i * nc + j].imag(), ref.imag());
        }
    }
}
