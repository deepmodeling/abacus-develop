#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/module_mixing/broyden_mixing.h"
#include "module_ri/Mix_DMk_2D.h"
#include "module_ri/Mix_Matrix.h"

/************************************************
 *  unit test of charge_mixing.cpp & Mix_DMk_2D.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Mix_Matrix::mix:
 *      mix the matrix data according to the set mixing mode.
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
        mix_data = std::vector<ModuleBase::matrix>(2);
        mix_complexdata = std::vector<ModuleBase::ComplexMatrix>(3);
        mix_data[0].create(nr, nc);
        mix_data[1].create(nr, nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_data[0](i, j) = i * nc + j;
                mix_data[1](i, j) = i * nc + j + 0.2;
            }
        }
        mix_complexdata[0].create(nr, nc);
        mix_complexdata[1].create(nr, nc);
        mix_complexdata[2].create(nr, nc);
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                mix_complexdata[0](i, j) = std::complex<double>{double(i), double(j)};
                mix_complexdata[1](i, j) = std::complex<double>{double(i), double(j) + 0.2};
                mix_complexdata[2](i, j) = std::complex<double>{double(i) + 0.8, double(j)};
            }
        }
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
    std::vector<ModuleBase::matrix> mix_data;
    std::vector<ModuleBase::ComplexMatrix> mix_complexdata;
    std::vector<std::vector<double>> mix_data_vector;
    std::vector<std::vector<std::complex<double>>> mix_complexdata_vector;
};

TEST_F(DM_Mixing_Test, Mix_Matrix)
{
    // Separete loop
    // matrix
    Mix_Matrix<ModuleBase::matrix> mix_matrix;
    mix_matrix.init(nullptr);
    mix_matrix.mixing_beta = mixing_beta;
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_matrix.mix(mix_data[istep], (istep == 0));
    }

    ModuleBase::matrix data_out = mix_matrix.get_data_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(data_out(i, j), (1 - mixing_beta) * mix_data[0](i, j) + mixing_beta * mix_data[1](i, j));
        }
    }
    // vector<double>
    Mix_Matrix<std::vector<double>> mix_vector;
    mix_vector.init(nullptr);
    mix_vector.mixing_beta = mixing_beta;
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_vector.mix(mix_data_vector[istep], (istep == 0));
    }

    std::vector<double> data_out_vector = mix_vector.get_data_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(data_out_vector[i * nc + j], (1 - mixing_beta) * mix_data_vector[0][i * nc + j] + mixing_beta * mix_data_vector[1][i * nc + j]);
        }
    }
    //vector<complex<double>>
    Mix_Matrix<std::vector<std::complex<double>>> mix_vector_complex;
    mix_vector_complex.init(nullptr);
    mix_vector_complex.mixing_beta = mixing_beta;
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_vector_complex.mix(mix_complexdata_vector[istep], (istep == 0));
    }

    std::vector<std::complex<double>> data_out_vector_complex = mix_vector_complex.get_data_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(data_out_vector_complex[i * nc + j].real(), ((1 - mixing_beta) * mix_complexdata_vector[0][i * nc + j] + mixing_beta * mix_complexdata_vector[1][i * nc + j]).real());
            EXPECT_DOUBLE_EQ(data_out_vector_complex[i * nc + j].imag(), ((1 - mixing_beta) * mix_complexdata_vector[0][i * nc + j] + mixing_beta * mix_complexdata_vector[1][i * nc + j]).imag());
        }
    }

    // Broyden mix
    Mix_Matrix<ModuleBase::ComplexMatrix> mix_complexmatrix;
    mix_complexmatrix.init(mixing);
    mixing->coef = { 1.1, -0.1 };
    for (int istep = 0; istep < 3; ++istep)
    {
        mix_complexmatrix.mix(mix_complexdata[istep], (istep == 0));
    }
    ModuleBase::ComplexMatrix com_data_out = mix_complexmatrix.get_data_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            std::complex<double> first_step_result
                = (1 - mixing_beta) * mix_complexdata[0](i, j) + mixing_beta * mix_complexdata[1](i, j);
            std::complex<double> second_step_result
                = (1 - mixing_beta) * first_step_result + mixing_beta * mix_complexdata[2](i, j);
            std::complex<double> ref = second_step_result * mixing->coef[1] + first_step_result * mixing->coef[0];
            EXPECT_DOUBLE_EQ(com_data_out(i, j).real(), ref.real());
            EXPECT_DOUBLE_EQ(com_data_out(i, j).imag(), ref.imag());
        }
    }
}

TEST_F(DM_Mixing_Test, Mix_DMk_2D)
{
    //Gamma only
    Mix_DMk_2D mix_dmk_gamma;
    mix_dmk_gamma.set_nks(1, true);
    mix_dmk_gamma.set_mixing(nullptr);
    mix_dmk_gamma.set_mixing_beta(mixing_beta);
    std::vector<std::vector<std::vector<double>>> dm_gamma(2);
    dm_gamma[0] = std::vector<std::vector<double>>(1);
    dm_gamma[0][0] = mix_data_vector[0];
    dm_gamma[1] = std::vector<std::vector<double>>(1);
    dm_gamma[1][0] = mix_data_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk_gamma.mix(dm_gamma[istep], (istep == 0));
    }
    std::vector<const std::vector<double>*> dm_gamma_out = mix_dmk_gamma.get_DMk_gamma_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            EXPECT_DOUBLE_EQ(dm_gamma_out[0][0][i * nc + j], (1 - mixing_beta) * mix_data_vector[0][i * nc + j] + mixing_beta * mix_data_vector[1][i * nc + j]);
        }
    }

    // not Gamma only
    Mix_DMk_2D mix_dmk;
    mix_dmk.set_nks(1, false);
    mix_dmk.set_mixing(nullptr);
    mix_dmk.set_mixing_beta(mixing_beta);
    std::vector<std::vector<std::vector<std::complex<double>>>> dm(2);
    dm[0] = std::vector<std::vector<std::complex<double>>>(1);
    dm[0][0] = mix_complexdata_vector[0];
    dm[1] = std::vector<std::vector<std::complex<double>>>(1);
    dm[1][0] = mix_complexdata_vector[1];
    for (int istep = 0; istep < 2; ++istep)
    {
        mix_dmk.mix(dm[istep], (istep == 0));
    }
    std::vector<const std::vector<std::complex<double>>*> dm_out = mix_dmk.get_DMk_k_out();
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            std::complex<double> first_step_result
                = (1 - mixing_beta) * mix_complexdata_vector[0][i * nc + j] + mixing_beta * mix_complexdata_vector[1][i * nc + j];
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].real(), first_step_result.real());
            EXPECT_DOUBLE_EQ(dm_out[0][0][i * nc + j].imag(), first_step_result.imag());
        }
    }
}