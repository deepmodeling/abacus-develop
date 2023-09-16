#include"gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_ri/Mix_Matrix.h"
#include "module_ri/Mix_DMk_2D.h"

class Mock_Charge_Mixing
{
public:
    MOCK_METHOD0(get_mixing_ndim, int());
};

enum class Mixing_Mode{No, Plain, Pulay};

class Charge_Mixing_Helper
{
public:
    Charge_Mixing_Helper(double mixing_beta, int mixing_ndim) :
        mixing_beta(mixing_beta), mixing_ndim(mixing_ndim){}
    ~Charge_Mixing_Helper() {}
    double mixing_beta;
    int mixing_ndim;
    int get_mixing_ndim() const { return this->mixing_ndim; }
};

class DM_Mixing_Test : public testing::Test
{
protected:
    int nbasis = 2;
    int maxiter = 6;
    // params for pulay mixing
    double mixing_beta = 0.4;
    int mixing_ndim = 4;
    int dstep = 3;
    double alpha[3] = {0.25, -0.5, 1 };

    void init_test();

    std::map<Mixing_Mode, std::vector<ModuleBase::matrix>> data_ref_double
        = {std::make_pair(Mixing_Mode::No, std::vector<ModuleBase::matrix>{ ModuleBase::matrix(nbasis, nbasis)}),
        std::make_pair(Mixing_Mode::Plain, std::vector<ModuleBase::matrix>{ ModuleBase::matrix(nbasis, nbasis)}),
        std::make_pair(Mixing_Mode::Pulay, std::vector<ModuleBase::matrix>{ ModuleBase::matrix(nbasis, nbasis)})};
     
    std::map<Mixing_Mode, std::vector<ModuleBase::ComplexMatrix>> data_ref_complex
        = { std::make_pair(Mixing_Mode::No, std::vector<ModuleBase::ComplexMatrix>{ ModuleBase::ComplexMatrix(nbasis, nbasis)}),
        std::make_pair(Mixing_Mode::Plain, std::vector<ModuleBase::ComplexMatrix>{ ModuleBase::ComplexMatrix(nbasis, nbasis)}),
        std::make_pair(Mixing_Mode::Pulay, std::vector<ModuleBase::ComplexMatrix>{ ModuleBase::ComplexMatrix(nbasis, nbasis)}) };

    // for distinguishing different iters
    ModuleBase::matrix mat1_double = ModuleBase::matrix(nbasis, nbasis);
    ModuleBase::ComplexMatrix mat1_complex = ModuleBase::ComplexMatrix(nbasis, nbasis);

    // for class Mix_Matrix
    Mix_Matrix<ModuleBase::matrix> mix_data_double;
    Mix_Matrix<ModuleBase::ComplexMatrix> mix_data_complex;

    // for pulay mixing
    std::vector<ModuleBase::matrix> simplemix_double;
    std::vector<ModuleBase::ComplexMatrix> simplemix_complex;
};

void DM_Mixing_Test::init_test()
{
    // set mixing_beta
    this->mix_data_double.mixing_beta = this->mixing_beta;
    this->mix_data_complex.mixing_beta = this->mixing_beta;
    // allocate memory for each iter
    for (auto mode : { Mixing_Mode::No, Mixing_Mode::Plain, Mixing_Mode::Pulay })
    {
        this->data_ref_double[mode].resize(this->maxiter);
        this->data_ref_double[mode][0] = ModuleBase::matrix(nbasis, nbasis);
        this->data_ref_complex[mode].resize(this->maxiter);
        this->data_ref_complex[mode][0] = ModuleBase::ComplexMatrix(nbasis, nbasis);
        for (int i = 0;i < nbasis;++i)
        {
            for(int j=0;j<nbasis;++j)
            {
                this->data_ref_double[mode][0](i, j) = 0.1 * (i * nbasis + j + 1);
                this->data_ref_complex[mode][0](i, j) = std::complex<double>(0.1 * (i * nbasis + j + 1), -0.1 * (i * nbasis + j + 1));
                this->mat1_double(i, j) = 0.01;
                this->mat1_complex(i, j) = std::complex<double>(0.01, -0.01);
            }
        }
    }
    
    // init data_ref
    for (int iter = 1;iter < this->maxiter;++iter)
    {
        ModuleBase::matrix& dm_in_double = this->data_ref_double[Mixing_Mode::No][iter];
        ModuleBase::ComplexMatrix& dm_in_complex = this->data_ref_complex[Mixing_Mode::No][iter];

        // No: dm_in for each iter
        dm_in_double = this->data_ref_double[Mixing_Mode::No][iter - 1] + mat1_double;
        dm_in_complex = this->data_ref_complex[Mixing_Mode::No][iter - 1] + mat1_complex;
        // Plain: beta*dm_in + (1-beta)*dm_out for each iter
        this->data_ref_double[Mixing_Mode::Plain][iter] = dm_in_double * this->mixing_beta
            + (1 - this->mixing_beta) * this->data_ref_double[Mixing_Mode::Plain][iter-1];
        this->data_ref_complex[Mixing_Mode::Plain][iter] = dm_in_complex * this->mixing_beta
            + (1 - this->mixing_beta) * this->data_ref_complex[Mixing_Mode::Plain][iter-1];
        //Pulay: directly using alpha
        this->simplemix_double.push_back(dm_in_double * this->mixing_beta
            + (1 - this->mixing_beta) * this->data_ref_double[Mixing_Mode::Pulay][iter - 1]);
        this->simplemix_complex.push_back(dm_in_complex * this->mixing_beta
            + (1 - this->mixing_beta) * this->data_ref_complex[Mixing_Mode::Pulay][iter - 1]);
        this->data_ref_double[Mixing_Mode::Pulay][iter] = this->simplemix_double[iter - 1];
        this->data_ref_complex[Mixing_Mode::Pulay][iter] = this->simplemix_complex[iter - 1];
        int mix_size = std::min(iter - 1, this->mixing_ndim);
        int tmp_idstep = (iter + 1) % dstep;
        for (int m = iter - 2;m > iter - 2 - mix_size;--m)
        {
            int alpha_index = (tmp_idstep + (dstep - mix_size)+ m ) % dstep;  // dstep = mixing_ndim -1
            this->data_ref_double[Mixing_Mode::Pulay][iter] += this->alpha[alpha_index] * (this->simplemix_double[m + 1] - this->simplemix_double[m]);
            this->data_ref_complex[Mixing_Mode::Pulay][iter] += this->alpha[alpha_index] * (this->simplemix_complex[m + 1] - this->simplemix_complex[m]);
        }
    }
}


// for class Mixing_data
TEST_F(DM_Mixing_Test, MixDataTest)
{
    this->init_test();

    // init mock
    Charge_Mixing_Helper chr_mix(mixing_beta, mixing_ndim);
    Mock_Charge_Mixing mock_chr_mix;

    EXPECT_CALL(mock_chr_mix, get_mixing_ndim()).WillRepeatedly(testing::Invoke(&chr_mix, &Charge_Mixing_Helper::get_mixing_ndim));

    //1. double
    for (auto mode : { Mixing_Mode::No, Mixing_Mode::Plain, Mixing_Mode::Pulay })
    {
        for (int iter = 1;iter < this->maxiter;++iter)
        {
            bool flag_restart = (iter == 1);
            this->mix_data_double.mix(data_ref_double[Mixing_Mode::No][iter - 1], flag_restart);
            ModuleBase::matrix data_out_double = mix_data_double.get_data_out();
            for (int i = 0;i < nbasis;++i)
                for (int j = 0;j < nbasis;++j)
                    EXPECT_NEAR(data_out_double(i, j), data_ref_double[mode][iter - 1](i, j), 1e-8);
        }
    }
    // 2. complex
    for (auto mode : { Mixing_Mode::No, Mixing_Mode::Plain, Mixing_Mode::Pulay })
    {
        for (int iter = 1;iter < this->maxiter;++iter)
        {
            bool flag_restart = (iter == 1);
            this->mix_data_complex.mix(data_ref_complex[Mixing_Mode::No][iter - 1], flag_restart);
            ModuleBase::ComplexMatrix data_out_complex = mix_data_complex.get_data_out();
            for (int i = 0;i < nbasis;++i)
                for (int j = 0;j < nbasis;++j)
                {
                    EXPECT_NEAR(data_out_complex(i, j).real(), data_ref_complex[mode][iter - 1](i, j).real(), 1e-8);
                    EXPECT_NEAR(data_out_complex(i, j).imag(), data_ref_complex[mode][iter - 1](i, j).imag(), 1e-8);
                }
                    
        }
    }
    
}

TEST_F(DM_Mixing_Test, MixDMk2D)
{
    this->init_test();
    
    // init mock
    Charge_Mixing_Helper chr_mix(mixing_beta, mixing_ndim);
    Mock_Charge_Mixing mock_chr_mix;

    EXPECT_CALL(mock_chr_mix, get_mixing_ndim()).WillRepeatedly(testing::Invoke(&chr_mix, &Charge_Mixing_Helper::get_mixing_ndim));
    
    //1. double
    Mix_DMk_2D mix_dm_gamma_2d;
    int nks = 1;
    mix_dm_gamma_2d.set_nks(nks, true);

    for (auto mode : { Mixing_Mode::No, Mixing_Mode::Plain, Mixing_Mode::Pulay })
    {
        mix_dm_gamma_2d.set_mixing_beta(this->mixing_beta);
        for (int iter = 1;iter < this->maxiter;++iter)
        {

            std::vector<ModuleBase::matrix> dm_gamma_in(nks);
            for (int ik = 0;ik < nks;++ik) dm_gamma_in[ik] = this->data_ref_double[Mixing_Mode::No][iter-1];

            bool flag_restart = (iter == 1);
            mix_dm_gamma_2d.mix(dm_gamma_in, flag_restart);
            
            std::vector<const ModuleBase::matrix*> dm_gamma_out = mix_dm_gamma_2d.get_DMk_gamma_out();
            for (int ik = 0;ik < nks;++ik)
                for (int i = 0;i < nbasis;++i)
                    for (int j = 0;j < nbasis;++j)
                        EXPECT_NEAR((*dm_gamma_out[ik])(i, j), data_ref_double[mode][iter - 1](i, j), 1e-8);
        }
    }

    //2. complex
    Mix_DMk_2D mix_dm_k_2d;
    nks = 2;
    mix_dm_k_2d.set_nks(nks, false);

    for (auto mode : { Mixing_Mode::No, Mixing_Mode::Plain, Mixing_Mode::Pulay })
    {
        mix_dm_k_2d.set_mixing_beta(this->mixing_beta);
        for (int iter = 1;iter < this->maxiter;++iter)
        {

            std::vector<ModuleBase::ComplexMatrix> dm_k_in(nks);
            for (int ik = 0;ik < nks;++ik) dm_k_in[ik] = this->data_ref_complex[Mixing_Mode::No][iter - 1];

            bool flag_restart = (iter == 1);
            mix_dm_k_2d.mix(dm_k_in, flag_restart);

            std::vector<const ModuleBase::ComplexMatrix*> dm_k_out = mix_dm_k_2d.get_DMk_k_out();
            for (int ik = 0;ik < nks;++ik)
                for (int i = 0;i < nbasis;++i)
                    for (int j = 0;j < nbasis;++j)
                    {
                        EXPECT_NEAR((*dm_k_out[ik])(i, j).real(), data_ref_complex[mode][iter - 1](i, j).real(), 1e-8);
                        EXPECT_NEAR((*dm_k_out[ik])(i, j).imag(), data_ref_complex[mode][iter - 1](i, j).imag(), 1e-8);
                    }
        }
    }

}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}