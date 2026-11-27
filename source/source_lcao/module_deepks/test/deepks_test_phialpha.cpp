#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_check.h"
#include "source_lcao/module_deepks/deepks_phialpha.h"

#include <complex>

template <typename T>
void test_deepks<T>::check_phialpha()
{
    std::vector<int> na;
    na.resize(ucell.ntype);
    for (int it = 0; it < ucell.ntype; it++)
    {
        na[it] = ucell.atoms[it].na;
    }
    this->ld.init(ORB, ucell.nat, ucell.ntype, kv.get_nkstot(), ParaO, na, GlobalV::ofs_running);

    DeePKS_domain::allocate_phialpha(this->cal_force, ucell, ORB, Test_Deepks::GridD, &ParaO, this->ld.phialpha);
    DeePKS_domain::build_phialpha(this->cal_force,
                                  ucell,
                                  ORB,
                                  Test_Deepks::GridD,
                                  &ParaO,
                                  overlap_orb_alpha_,
                                  this->ld.phialpha);
    DeePKS_domain::check_phialpha(this->cal_force,
                                  ucell,
                                  ORB,
                                  Test_Deepks::GridD,
                                  &ParaO,
                                  this->ld.phialpha,
                                  0);

    this->assert_file_matches_reference("phialpha.dat", "phialpha_ref.dat");
    this->assert_file_matches_reference("dphialpha_x.dat", "dphialpha_x_ref.dat");
    this->assert_file_matches_reference("dphialpha_y.dat", "dphialpha_y_ref.dat");
    this->assert_file_matches_reference("dphialpha_z.dat", "dphialpha_z_ref.dat");
}

template void test_deepks<double>::check_phialpha();
template void test_deepks<std::complex<double>>::check_phialpha();

template <typename T>
void run_deepks_unit_phialpha(test_deepks<T>& test)
{
    test.check_phialpha();
}

template void run_deepks_unit_phialpha<double>(test_deepks<double>& test);
template void run_deepks_unit_phialpha<std::complex<double>>(test_deepks<std::complex<double>>& test);
