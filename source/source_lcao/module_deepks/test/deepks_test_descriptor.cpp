#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_descriptor.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_descriptor(std::vector<torch::Tensor>& descriptor)
{
    DeePKS_domain::cal_descriptor(ucell.nat, this->ld.deepks_param, this->ld.pdm, descriptor);
    DeePKS_domain::check_descriptor(this->ld.deepks_param, ucell, "./", descriptor, 0);
    this->assert_file_matches_reference("deepks_desc.dat", "descriptor_ref.dat");
}

template void test_deepks<double>::check_descriptor(std::vector<torch::Tensor>& descriptor);
template void test_deepks<std::complex<double>>::check_descriptor(std::vector<torch::Tensor>& descriptor);

template <typename T>
void run_deepks_unit_descriptor(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
}

template void run_deepks_unit_descriptor<double>(test_deepks<double>& test);
template void run_deepks_unit_descriptor<std::complex<double>>(test_deepks<std::complex<double>>& test);
