#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_orbital.h"

#include <complex>
#include <iomanip>

template <typename T>
void test_deepks<T>::check_o_delta()
{
    const int nks = kv.get_nkstot();
    ModuleBase::matrix o_delta;
    o_delta.create(nks, 1);
    DeePKS_domain::cal_o_delta<T>(dm, ld.V_delta, o_delta, ParaO, nks, this->nspin);
    std::ofstream ofs("o_delta.dat");
    ofs << std::setprecision(10);
    o_delta.print(ofs);
    ofs.close();
    this->assert_file_matches_reference("o_delta.dat", "o_delta_ref.dat");
}

template void test_deepks<double>::check_o_delta();
template void test_deepks<std::complex<double>>::check_o_delta();

template <typename T>
void run_deepks_unit_o_delta(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_edelta(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_e_deltabands();
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_o_delta();
}

template void run_deepks_unit_o_delta<double>(test_deepks<double>& test);
template void run_deepks_unit_o_delta<std::complex<double>>(test_deepks<std::complex<double>>& test);
