#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_basic.h"
#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_vdpre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_vdpre()
{
    std::vector<torch::Tensor> gevdm;
    torch::Tensor vdpre;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.deepks_param, this->ld.pdm, gevdm);
    DeePKS_domain::cal_v_delta_precalc<T>(this->nlocal,
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
                                          vdpre);
    DeePKS_domain::check_tensor<T>(vdpre, "v_delta_precalc.dat", 0);
    this->assert_file_matches_reference("v_delta_precalc.dat", "vdpre_ref.dat");
}

template void test_deepks<double>::check_vdpre();
template void test_deepks<std::complex<double>>::check_vdpre();

template <typename T>
void run_deepks_unit_vdpre(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_vdpre();
}

template void run_deepks_unit_vdpre<double>(test_deepks<double>& test);
template void run_deepks_unit_vdpre<std::complex<double>>(test_deepks<std::complex<double>>& test);
