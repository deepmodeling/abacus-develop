#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_force.h"

#include <complex>
#include <iomanip>

template <typename T>
void test_deepks<T>::check_f_delta_and_stress_delta()
{
    ModuleBase::matrix fvnl_dalpha;
    fvnl_dalpha.create(ucell.nat, 3);

    ModuleBase::matrix svnl_dalpha;
    svnl_dalpha.create(3, 3);
    const int cal_stress = 1;
    const int nks = kv.get_nkstot();
    DeePKS_domain::cal_f_delta<T>(ucell,
                                  ORB,
                                  Test_Deepks::GridD,
                                  ParaO,
                                  nks,
                                  this->ld.deepks_param,
                                  kv.kvec_d,
                                  this->ld.phialpha,
                                  fvnl_dalpha,
                                  cal_stress,
                                  svnl_dalpha,
                                  this->ld.dm_r,
                                  this->ld.gedm);
    std::ofstream ofs_f("F_delta.dat");
    std::ofstream ofs_s("stress_delta.dat");
    ofs_f << std::setprecision(10);
    ofs_s << std::setprecision(10);
    fvnl_dalpha.print(ofs_f);
    ofs_f.close();
    svnl_dalpha.print(ofs_s);
    ofs_s.close();

    this->assert_file_matches_reference("F_delta.dat", "F_delta_ref.dat");
    this->assert_file_matches_reference("stress_delta.dat", "stress_delta_ref.dat");
}

template void test_deepks<double>::check_f_delta_and_stress_delta();
template void test_deepks<std::complex<double>>::check_f_delta_and_stress_delta();

template <typename T>
void run_deepks_unit_force_stress_delta(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_edelta(test, descriptor);
    if (testing::Test::HasFatalFailure())
    {
        return;
    }
    test.check_f_delta_and_stress_delta();
}

template void run_deepks_unit_force_stress_delta<double>(test_deepks<double>& test);
template void run_deepks_unit_force_stress_delta<std::complex<double>>(test_deepks<std::complex<double>>& test);
