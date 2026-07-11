#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_pdm.h"

#include <complex>
#include <gtest/gtest.h>

template <typename T>
void test_deepks<T>::read_dm(const int nks)
{
    dm.resize(nks);
    std::stringstream ss;
    for (int ik = 0; ik < nks; ik++)
    {
        ss.str("");
        if (nks == 1)
        {
            ss << "dm";
        }
        else
        {
            ss << "dm_" << ik;
        }
        std::ifstream ifs(ss.str().c_str());
        ASSERT_TRUE(ifs.is_open()) << "Cannot open density matrix file " << ss.str();
        dm[ik].create(this->nlocal, this->nlocal);

        for (int mu = 0; mu < this->nlocal; mu++)
        {
            for (int nu = 0; nu < this->nlocal; nu++)
            {
                T c;
                ASSERT_TRUE(ifs >> c) << "Failed to read " << ss.str() << " at (" << mu << ", " << nu << ")";
                dm[ik](mu, nu) = c;
            }
        }
    }
}

template <typename T>
void test_deepks<T>::set_dm_new()
{
    dm_new.resize(dm.size());
    for (int i = 0; i < dm.size(); i++)
    {
        dm_new[i].resize(dm[i].nr * dm[i].nc);
        dm_new[i].assign(dm[i].c, dm[i].c + dm[i].nr * dm[i].nc);
    }
}

template <typename T>
void test_deepks<T>::set_p_elec_DM()
{
    int nk = 1;
    if (this->gamma_only_local)
    {
        nk = this->nspin;
        this->p_elec_DM = new elecstate::DensityMatrix<T, double>(&ParaO, this->nspin);
    }
    else
    {
        nk = kv.get_nkstot();
        this->p_elec_DM
            = new elecstate::DensityMatrix<T, double>(&ParaO, this->nspin, kv.kvec_d, kv.get_nkstot() / this->nspin);
    }
    p_elec_DM->init_DMR(&Test_Deepks::GridD, &ucell);

    for (int ik = 0; ik < nk; ik++)
    {
        p_elec_DM->set_DMK_pointer(ik, dm_new[ik].data());
    }
    p_elec_DM->cal_DMR();
}

template <typename T>
void test_deepks<T>::check_pdm()
{
    this->read_dm(kv.get_nkstot());
    this->set_dm_new();
    this->set_p_elec_DM();
    this->ld.init_DMR(ucell, ORB, ParaO, Test_Deepks::GridD);
    DeePKS_domain::update_dmr(kv.kvec_d,
                              p_elec_DM->get_DMK_vector(),
                              ucell,
                              ORB,
                              ParaO,
                              Test_Deepks::GridD,
                              this->ld.dm_r);
    DeePKS_domain::cal_pdm<T>(this->ld.init_pdm,
                              this->ld.deepks_param,
                              kv.kvec_d,
                              this->ld.dm_r,
                              this->ld.phialpha,
                              ucell,
                              ORB,
                              Test_Deepks::GridD,
                              ParaO,
                              this->ld.pdm);
    DeePKS_domain::check_pdm(this->ld.deepks_param, this->ld.pdm);
    this->assert_file_matches_reference("deepks_projdm.dat", "pdm_ref.dat");
}

template void test_deepks<double>::read_dm(const int nks);
template void test_deepks<std::complex<double>>::read_dm(const int nks);
template void test_deepks<double>::set_dm_new();
template void test_deepks<std::complex<double>>::set_dm_new();
template void test_deepks<double>::set_p_elec_DM();
template void test_deepks<std::complex<double>>::set_p_elec_DM();
template void test_deepks<double>::check_pdm();
template void test_deepks<std::complex<double>>::check_pdm();

template <typename T>
void run_deepks_unit_pdm(test_deepks<T>& test)
{
    DeepksTestRunner::build_pdm(test);
}

template void run_deepks_unit_pdm<double>(test_deepks<double>& test);
template void run_deepks_unit_pdm<std::complex<double>>(test_deepks<std::complex<double>>& test);
