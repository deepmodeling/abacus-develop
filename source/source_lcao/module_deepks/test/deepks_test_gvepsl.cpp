#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_basic.h"
#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_spre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_gvepsl(torch::Tensor& gdmepsl)
{
    std::vector<torch::Tensor> gevdm;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.deepks_param, this->ld.pdm, gevdm);
    torch::Tensor gvepsl;
    DeePKS_domain::cal_gvepsl(ucell.nat, this->ld.deepks_param, gevdm, gdmepsl, gvepsl, 0);
    DeePKS_domain::check_tensor<double>(gvepsl, "gvepsl.dat", 0);
    this->assert_file_matches_reference("gvepsl.dat", "gvepsl_ref.dat");
}

template void test_deepks<double>::check_gvepsl(torch::Tensor& gdmepsl);
template void test_deepks<std::complex<double>>::check_gvepsl(torch::Tensor& gdmepsl);

template <typename T>
void run_deepks_unit_gvepsl(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    torch::Tensor gdmepsl;
    test.check_gdmepsl(gdmepsl);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_gvepsl(gdmepsl);
}

template void run_deepks_unit_gvepsl<double>(test_deepks<double>& test);
template void run_deepks_unit_gvepsl<std::complex<double>>(test_deepks<std::complex<double>>& test);
