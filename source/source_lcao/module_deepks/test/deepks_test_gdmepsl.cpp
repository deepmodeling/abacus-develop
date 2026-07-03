#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_spre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_gdmepsl(torch::Tensor& gdmepsl)
{
    DeePKS_domain::cal_gdmepsl<T>(kv.get_nkstot(),
                                  this->ld.deepks_param,
                                  kv.kvec_d,
                                  this->ld.phialpha,
                                  this->ld.dm_r,
                                  ucell,
                                  ORB,
                                  ParaO,
                                  Test_Deepks::GridD,
                                  gdmepsl);
    DeePKS_domain::check_tensor<double>(gdmepsl, "gdmepsl.dat", 0);
    this->assert_file_matches_reference("gdmepsl.dat", "gdmepsl_ref.dat");
}

template void test_deepks<double>::check_gdmepsl(torch::Tensor& gdmepsl);
template void test_deepks<std::complex<double>>::check_gdmepsl(torch::Tensor& gdmepsl);

template <typename T>
void run_deepks_unit_gdmepsl(test_deepks<T>& test)
{
    DeepksTestRunner::build_pdm(test);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    torch::Tensor gdmepsl;
    test.check_gdmepsl(gdmepsl);
}

template void run_deepks_unit_gdmepsl<double>(test_deepks<double>& test);
template void run_deepks_unit_gdmepsl<std::complex<double>>(test_deepks<std::complex<double>>& test);
