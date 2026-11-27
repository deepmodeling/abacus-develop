#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_basic.h"
#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_orbpre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_orbpre()
{
    using TH = std::conditional_t<std::is_same<T, double>::value, ModuleBase::matrix, ModuleBase::ComplexMatrix>;
    std::vector<torch::Tensor> gevdm;
    torch::Tensor orbpre;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.deepks_param, this->ld.pdm, gevdm);
    DeePKS_domain::cal_orbital_precalc<T, TH>(dm,
                                              ucell.nat,
                                              kv.get_nkstot(),
                                              this->ld.deepks_param,
                                              kv.kvec_d,
                                              this->ld.phialpha,
                                              gevdm,
                                              ucell,
                                              ORB,
                                              ParaO,
                                              Test_Deepks::GridD,
                                              orbpre);
    DeePKS_domain::check_tensor<double>(orbpre, "orbital_precalc.dat", 0);
    this->assert_file_matches_reference("orbital_precalc.dat", "orbpre_ref.dat");
}

template void test_deepks<double>::check_orbpre();
template void test_deepks<std::complex<double>>::check_orbpre();

template <typename T>
void run_deepks_unit_orbpre(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_orbpre();
}

template void run_deepks_unit_orbpre<double>(test_deepks<double>& test);
template void run_deepks_unit_orbpre<std::complex<double>>(test_deepks<std::complex<double>>& test);
