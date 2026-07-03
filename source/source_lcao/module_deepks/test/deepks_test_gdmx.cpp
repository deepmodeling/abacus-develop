#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_fpre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_gdmx(torch::Tensor& gdmx)
{
    DeePKS_domain::cal_gdmx<T>(kv.get_nkstot(),
                               this->ld.deepks_param,
                               kv.kvec_d,
                               this->ld.phialpha,
                               this->ld.dm_r,
                               ucell,
                               ORB,
                               ParaO,
                               Test_Deepks::GridD,
                               gdmx);
    DeePKS_domain::check_tensor<double>(gdmx, "gdmx.dat", 0);
    this->assert_file_matches_reference("gdmx.dat", "gdmx_ref.dat");
}

template void test_deepks<double>::check_gdmx(torch::Tensor& gdmx);
template void test_deepks<std::complex<double>>::check_gdmx(torch::Tensor& gdmx);

template <typename T>
void run_deepks_unit_gdmx(test_deepks<T>& test)
{
    DeepksTestRunner::build_pdm(test);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    torch::Tensor gdmx;
    test.check_gdmx(gdmx);
}

template void run_deepks_unit_gdmx<double>(test_deepks<double>& test);
template void run_deepks_unit_gdmx<std::complex<double>>(test_deepks<std::complex<double>>& test);
