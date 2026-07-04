#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_basic.h"
#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_fpre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_gvx(torch::Tensor& gdmx)
{
    std::vector<torch::Tensor> gevdm;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.deepks_param, this->ld.pdm, gevdm);
    torch::Tensor gvx;
    DeePKS_domain::cal_gvx(ucell.nat, this->ld.deepks_param, gevdm, gdmx, gvx, 0);
    DeePKS_domain::check_tensor<double>(gvx, "gvx.dat", 0);
    this->assert_file_matches_reference("gvx.dat", "gvx_ref.dat");
}

template void test_deepks<double>::check_gvx(torch::Tensor& gdmx);
template void test_deepks<std::complex<double>>::check_gvx(torch::Tensor& gdmx);

template <typename T>
void run_deepks_unit_gvx(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    torch::Tensor gdmx;
    test.check_gdmx(gdmx);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_gvx(gdmx);
}

template void run_deepks_unit_gvx<double>(test_deepks<double>& test);
template void run_deepks_unit_gvx<std::complex<double>>(test_deepks<std::complex<double>>& test);
