#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_basic.h"
#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_vdrpre.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_vdrpre()
{
    std::vector<torch::Tensor> gevdm;
    torch::Tensor vdrpre;
    torch::Tensor overlap_out;
    torch::Tensor iRmat;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.deepks_param, this->ld.pdm, gevdm);
    const int R_size = 3;
    DeePKS_domain::cal_vdr_precalc(this->nlocal,
                                   ucell.nat,
                                   kv.get_nkstot(),
                                   R_size,
                                   this->ld.deepks_param,
                                   kv.kvec_d,
                                   this->ld.phialpha,
                                   gevdm,
                                   ucell,
                                   ORB,
                                   ParaO,
                                   Test_Deepks::GridD,
                                   vdrpre);
    DeePKS_domain::prepare_phialpha_iRmat(this->nlocal,
                                          R_size,
                                          this->ld.deepks_param,
                                          this->ld.phialpha,
                                          ucell,
                                          ORB,
                                          ParaO,
                                          Test_Deepks::GridD,
                                          overlap_out,
                                          iRmat);
    torch::Tensor vdrpre_sliced = vdrpre.slice(0, 0, 2, 1).slice(1, 0, 1, 1).slice(2, 0, 1, 1);
    DeePKS_domain::check_tensor<double>(vdrpre_sliced, "vdr_precalc.dat", 0);
    DeePKS_domain::check_tensor<double>(overlap_out, "phialpha_r.dat", 0);
    DeePKS_domain::check_tensor<int>(iRmat, "iRmat.dat", 0);
    this->assert_file_matches_reference("vdr_precalc.dat", "vdrpre_ref.dat");
    this->assert_file_matches_reference("phialpha_r.dat", "phialpha_r_ref.dat");
    this->assert_file_matches_reference("iRmat.dat", "iRmat_ref.dat");
}

template void test_deepks<double>::check_vdrpre();
template void test_deepks<std::complex<double>>::check_vdrpre();

template <typename T>
void run_deepks_unit_vdrpre(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_descriptor(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_vdrpre();
}

template void run_deepks_unit_vdrpre<double>(test_deepks<double>& test);
template void run_deepks_unit_vdrpre<std::complex<double>>(test_deepks<std::complex<double>>& test);
