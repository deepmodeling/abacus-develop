#include "deepks_test_runner.h"

#include "source_lcao/hs_matrix_k.hpp"
#include "source_lcao/module_operator_lcao/deepks_lcao.h"

#include <complex>
#include <iomanip>

template <typename T>
void test_deepks<T>::cal_V_delta()
{
    hamilt::HS_Matrix_K<T> hsk(&ParaO);
    hamilt::HContainer<double> hR(ucell, &ParaO);
    hamilt::DeePKS<hamilt::OperatorLCAO<T, double>> op_deepks(&hsk,
                                                              kv.kvec_d,
                                                              &hR,
                                                              &ucell,
                                                              &Test_Deepks::GridD,
                                                              &overlap_orb_alpha_,
                                                              &ORB,
                                                              kv.get_nkstot(),
                                                              p_elec_DM,
                                                              &this->ld);
    for (int ik = 0; ik < kv.get_nkstot(); ++ik)
    {
        op_deepks.init(ik);
    }
}

template <typename T>
void test_deepks<T>::check_e_deltabands()
{
    this->cal_V_delta();
    this->ld.dpks_cal_e_delta_band(dm_new, kv.get_nkstot());

    std::ofstream ofs("E_delta_bands.dat");
    ofs << std::setprecision(10) << this->ld.e_delta_band << std::endl;
    ofs.close();
    this->assert_file_matches_reference("E_delta_bands.dat", "E_delta_bands_ref.dat");
}

template void test_deepks<double>::cal_V_delta();
template void test_deepks<std::complex<double>>::cal_V_delta();
template void test_deepks<double>::check_e_deltabands();
template void test_deepks<std::complex<double>>::check_e_deltabands();

template <typename T>
void run_deepks_unit_e_deltabands(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_edelta(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_e_deltabands();
}

template void run_deepks_unit_e_deltabands<double>(test_deepks<double>& test);
template void run_deepks_unit_e_deltabands<std::complex<double>>(test_deepks<std::complex<double>>& test);
