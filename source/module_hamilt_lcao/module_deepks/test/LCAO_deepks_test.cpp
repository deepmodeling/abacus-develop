#include "LCAO_deepks_test.h"
#define private public
#include "module_parameter/parameter.h"

#include <torch/script.h>
#include <torch/torch.h>
#undef private
#include "module_hamilt_lcao/hamilt_lcaodft/hs_matrix_k.hpp"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/deepks_lcao.h"
namespace Test_Deepks
{
Grid_Driver GridD(PARAM.input.test_deconstructor, PARAM.input.test_grid);
}

test_deepks::test_deepks()
{
}

test_deepks::~test_deepks()
{
}

void test_deepks::check_dstable()
{
    // OGT.talpha.print_Table_DSR(ORB);
    // this->compare_with_ref("S_I_mu_alpha.dat","S_I_mu_alpha_ref.dat");
}

void test_deepks::check_phialpha()
{
    std::vector<int> na;
    na.resize(ucell.ntype);
    for (int it = 0; it < ucell.ntype; it++)
    {
        na[it] = ucell.atoms[it].na;
    }
    this->ld.init(ORB, ucell.nat, ucell.ntype, kv.nkstot, ParaO, na);

    DeePKS_domain::allocate_phialpha(PARAM.input.cal_force, ucell, ORB, Test_Deepks::GridD, &ParaO, this->ld.phialpha);

    DeePKS_domain::build_phialpha(PARAM.input.cal_force,
                                  ucell,
                                  ORB,
                                  Test_Deepks::GridD,
                                  &ParaO,
                                  overlap_orb_alpha_,
                                  this->ld.phialpha);

    DeePKS_domain::check_phialpha(PARAM.input.cal_force, ucell, ORB, Test_Deepks::GridD, &ParaO, this->ld.phialpha);

    this->compare_with_ref("phialpha.dat", "phialpha_ref.dat");
    this->compare_with_ref("dphialpha_x.dat", "dphialpha_x_ref.dat");
    this->compare_with_ref("dphialpha_y.dat", "dphialpha_y_ref.dat");
    this->compare_with_ref("dphialpha_z.dat", "dphialpha_z_ref.dat");
}

void test_deepks::read_dm()
{
    std::ifstream ifs("dm");
    dm.resize(1);
    dm[0].create(PARAM.sys.nlocal, PARAM.sys.nlocal);

    for (int mu = 0; mu < PARAM.sys.nlocal; mu++)
    {
        for (int nu = 0; nu < PARAM.sys.nlocal; nu++)
        {
            double c;
            ifs >> c;
            dm[0](mu, nu) = c;
        }
    }
}

void test_deepks::read_dm_k(const int nks)
{
    dm_k.resize(nks);
    std::stringstream ss;
    for (int ik = 0; ik < nks; ik++)
    {
        ss.str("");
        ss << "dm_" << ik;
        std::ifstream ifs(ss.str().c_str());
        dm_k[ik].create(PARAM.sys.nlocal, PARAM.sys.nlocal);

        for (int mu = 0; mu < PARAM.sys.nlocal; mu++)
        {
            for (int nu = 0; nu < PARAM.sys.nlocal; nu++)
            {
                std::complex<double> c;
                ifs >> c;
                dm_k[ik](mu, nu) = c;
            }
        }
    }
}

void test_deepks::set_dm_new()
{
    // dm_gamma
    dm_new.resize(dm.size());
    for (int i = 0; i < dm.size(); i++)
    {
        dm_new[i].resize(dm[i].nr * dm[i].nc);
        dm_new[i].assign(dm[i].c, dm[i].c + dm[i].nr * dm[i].nc);
    }
}

void test_deepks::set_dm_k_new()
{
    // dm_k
    dm_k_new.resize(dm_k.size());
    for (int i = 0; i < dm_k.size(); i++)
    {
        dm_k_new[i].resize(dm_k[i].nr * dm_k[i].nc);
        dm_k_new[i].assign(dm_k[i].c, dm_k[i].c + dm_k[i].nr * dm_k[i].nc);
    }
}

void test_deepks::set_p_elec_DM()
{
    // gamma
    int nspin = PARAM.inp.nspin;
    this->p_elec_DM = new elecstate::DensityMatrix<double, double>(&ParaO, nspin);
    p_elec_DM->init_DMR(&Test_Deepks::GridD, &ucell);
    for (int ik = 0; ik < nspin; ik++)
    {
        p_elec_DM->set_DMK_pointer(ik, dm_new[ik].data());
    }
    p_elec_DM->cal_DMR();
}

void test_deepks::set_p_elec_DM_k()
{
    // multi k
    this->p_elec_DM_k = new elecstate::DensityMatrix<std::complex<double>, double>(&ParaO,
                                                                                   PARAM.inp.nspin,
                                                                                   kv.kvec_d,
                                                                                   kv.nkstot / PARAM.inp.nspin);
    p_elec_DM_k->init_DMR(&Test_Deepks::GridD, &ucell);
    for (int ik = 0; ik < kv.nkstot; ++ik)
    {
        p_elec_DM_k->set_DMK_pointer(ik, dm_k_new[ik].data());
    }
    p_elec_DM_k->cal_DMR();
}

void test_deepks::check_pdm()
{
    if (PARAM.sys.gamma_only_local)
    {
        this->read_dm();
        this->set_dm_new();
        this->set_p_elec_DM();
        DeePKS_domain::cal_pdm(this->ld.init_pdm,
                               this->ld.inlmax,
                               this->ld.lmaxd,
                               this->ld.inl_l,
                               this->ld.inl_index,
                               p_elec_DM,
                               this->ld.phialpha,
                               ucell,
                               ORB,
                               Test_Deepks::GridD,
                               ParaO,
                               this->ld.pdm);
    }
    else
    {
        this->read_dm_k(kv.nkstot);
        this->set_dm_k_new();
        this->set_p_elec_DM_k();
        DeePKS_domain::cal_pdm(this->ld.init_pdm,
                               this->ld.inlmax,
                               this->ld.lmaxd,
                               this->ld.inl_l,
                               this->ld.inl_index,
                               p_elec_DM_k,
                               this->ld.phialpha,
                               ucell,
                               ORB,
                               Test_Deepks::GridD,
                               ParaO,
                               this->ld.pdm);
    }
    DeePKS_domain::check_pdm(this->ld.inlmax, this->ld.inl_l, this->ld.pdm);
    this->compare_with_ref("pdm.dat", "pdm_ref.dat");
}

void test_deepks::check_gdmx(torch::Tensor& gdmx)
{
    if (PARAM.sys.gamma_only_local)
    {
        DeePKS_domain::cal_gdmx(this->ld.lmaxd,
                                this->ld.inlmax,
                                kv.nkstot,
                                kv.kvec_d,
                                this->ld.phialpha,
                                this->ld.inl_index,
                                dm_new,
                                ucell,
                                ORB,
                                ParaO,
                                Test_Deepks::GridD,
                                gdmx);
    }
    else
    {
        DeePKS_domain::cal_gdmx(this->ld.lmaxd,
                                this->ld.inlmax,
                                kv.nkstot,
                                kv.kvec_d,
                                this->ld.phialpha,
                                this->ld.inl_index,
                                dm_k_new,
                                ucell,
                                ORB,
                                ParaO,
                                Test_Deepks::GridD,
                                gdmx);
    }
    DeePKS_domain::check_gdmx(gdmx);

    for (int ia = 0; ia < ucell.nat; ia++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gdmx_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmx_" << ia << "_ref.dat";

        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gdmy_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmy_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gdmz_" << ia << ".dat";
        ss1.str("");
        ss1 << "gdmz_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_gdmepsl(torch::Tensor& gdmepsl)
{
    if (PARAM.sys.gamma_only_local)
    {
        DeePKS_domain::cal_gdmepsl(this->ld.lmaxd,
                                   this->ld.inlmax,
                                   kv.nkstot,
                                   kv.kvec_d,
                                   this->ld.phialpha,
                                   this->ld.inl_index,
                                   dm_new,
                                   ucell,
                                   ORB,
                                   ParaO,
                                   Test_Deepks::GridD,
                                   gdmepsl);
    }
    else
    {
        DeePKS_domain::cal_gdmepsl(this->ld.lmaxd,
                                   this->ld.inlmax,
                                   kv.nkstot,
                                   kv.kvec_d,
                                   this->ld.phialpha,
                                   this->ld.inl_index,
                                   dm_k_new,
                                   ucell,
                                   ORB,
                                   ParaO,
                                   Test_Deepks::GridD,
                                   gdmepsl);
    }
    DeePKS_domain::check_gdmepsl(gdmepsl);

    for (int i = 0; i < 6; i++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gdmepsl_" << i << ".dat";
        ss1.str("");
        ss1 << "gdmepsl_" << i << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_descriptor(std::vector<torch::Tensor>& descriptor)
{
    DeePKS_domain::cal_descriptor(ucell.nat,
                                  this->ld.inlmax,
                                  this->ld.inl_l,
                                  this->ld.pdm,
                                  descriptor,
                                  this->ld.des_per_atom);
    DeePKS_domain::check_descriptor(this->ld.inlmax, this->ld.des_per_atom, this->ld.inl_l, ucell, "./", descriptor);
    this->compare_with_ref("deepks_desc.dat", "descriptor_ref.dat");
}

void test_deepks::check_gvx(torch::Tensor& gdmx)
{
    std::vector<torch::Tensor> gevdm;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.inlmax, this->ld.inl_l, this->ld.pdm, gevdm);
    torch::Tensor gvx;
    DeePKS_domain::cal_gvx(ucell.nat, this->ld.inlmax, this->ld.des_per_atom, this->ld.inl_l, gevdm, gdmx, gvx);
    DeePKS_domain::check_gvx(gvx);

    for (int ia = 0; ia < ucell.nat; ia++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gvx_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvx_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gvy_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvy_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());

        ss.str("");
        ss << "gvz_" << ia << ".dat";
        ss1.str("");
        ss1 << "gvz_" << ia << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_gvepsl(torch::Tensor& gdmepsl)
{
    std::vector<torch::Tensor> gevdm;
    DeePKS_domain::cal_gevdm(ucell.nat, this->ld.inlmax, this->ld.inl_l, this->ld.pdm, gevdm);
    torch::Tensor gvepsl;
    DeePKS_domain::cal_gvepsl(ucell.nat,
                              this->ld.inlmax,
                              this->ld.des_per_atom,
                              this->ld.inl_l,
                              gevdm,
                              gdmepsl,
                              gvepsl);
    DeePKS_domain::check_gvepsl(gvepsl);

    for (int i = 0; i < 6; i++)
    {
        std::stringstream ss;
        std::stringstream ss1;
        ss.str("");
        ss << "gvepsl_" << i << ".dat";
        ss1.str("");
        ss1 << "gvepsl_" << i << "_ref.dat";
        this->compare_with_ref(ss.str(), ss1.str());
    }
}

void test_deepks::check_edelta(std::vector<torch::Tensor>& descriptor)
{
    DeePKS_domain::load_model("model.ptg", ld.model_deepks);
    if (PARAM.sys.gamma_only_local)
    {
        ld.allocate_V_delta(ucell.nat, 1); // 1 for gamma-only
    }
    else
    {
        ld.allocate_V_delta(ucell.nat, kv.nkstot);
    }
    DeePKS_domain::cal_gedm(ucell.nat,
                            this->ld.lmaxd,
                            this->ld.nmaxd,
                            this->ld.inlmax,
                            this->ld.des_per_atom,
                            this->ld.inl_l,
                            descriptor,
                            this->ld.pdm,
                            this->ld.model_deepks,
                            this->ld.gedm,
                            this->ld.E_delta);

    std::ofstream ofs("E_delta.dat");
    ofs << std::setprecision(10) << this->ld.E_delta << std::endl;
    ofs.close();
    this->compare_with_ref("E_delta.dat", "E_delta_ref.dat");

    DeePKS_domain::check_gedm(this->ld.inlmax, this->ld.inl_l, this->ld.gedm);
    this->compare_with_ref("gedm.dat", "gedm_ref.dat");
}

void test_deepks::cal_H_V_delta()
{
    hamilt::HS_Matrix_K<double>* hsk = new hamilt::HS_Matrix_K<double>(&ParaO);
    hamilt::HContainer<double>* hR = new hamilt::HContainer<double>(ucell, &ParaO);
    hamilt::Operator<double>* op_deepks
        = new hamilt::DeePKS<hamilt::OperatorLCAO<double, double>>(hsk,
                                                                   kv.kvec_d,
                                                                   hR, // no explicit call yet
                                                                   &ucell,
                                                                   &Test_Deepks::GridD,
                                                                   &overlap_orb_alpha_,
                                                                   &ORB,
                                                                   kv.nkstot,
                                                                   p_elec_DM,
                                                                   &this->ld);
    for (int ik = 0; ik < kv.nkstot; ++ik)
    {
        op_deepks->init(ik);
    }
}

void test_deepks::cal_H_V_delta_k()
{
    hamilt::HS_Matrix_K<std::complex<double>>* hsk = new hamilt::HS_Matrix_K<std::complex<double>>(&ParaO);
    hamilt::HContainer<double>* hR = new hamilt::HContainer<double>(&ParaO);

    hamilt::Operator<std::complex<double>>* op_deepks
        = new hamilt::DeePKS<hamilt::OperatorLCAO<std::complex<double>, double>>(hsk,
                                                                                 kv.kvec_d,
                                                                                 hR, // no explicit call yet
                                                                                 &ucell,
                                                                                 &Test_Deepks::GridD,
                                                                                 &overlap_orb_alpha_,
                                                                                 &ORB,
                                                                                 kv.nkstot,
                                                                                 p_elec_DM_k,
                                                                                 &this->ld);
    for (int ik = 0; ik < kv.nkstot; ++ik)
    {
        op_deepks->init(ik);
    }
}

void test_deepks::check_e_deltabands()
{
    if (PARAM.sys.gamma_only_local)
    {
        this->cal_H_V_delta();
        this->ld.dpks_cal_e_delta_band(dm_new, 1);
    }
    else
    {
        this->cal_H_V_delta_k();
        this->ld.dpks_cal_e_delta_band(dm_k_new, kv.nkstot);
    }

    std::ofstream ofs("E_delta_bands.dat");
    ofs << std::setprecision(10) << this->ld.e_delta_band << std::endl;
    ofs.close();
    this->compare_with_ref("E_delta_bands.dat", "E_delta_bands_ref.dat");
}

void test_deepks::check_f_delta_and_stress_delta()
{
    ModuleBase::matrix fvnl_dalpha;
    fvnl_dalpha.create(ucell.nat, 3);

    ModuleBase::matrix svnl_dalpha;
    svnl_dalpha.create(3, 3);
    const int cal_stress = 1;
    if (PARAM.sys.gamma_only_local)
    {
        const int nks = 1;
        DeePKS_domain::cal_f_delta<double>(dm_new,
                                           ucell,
                                           ORB,
                                           Test_Deepks::GridD,
                                           ParaO,
                                           nks,
                                           kv.kvec_d,
                                           this->ld.phialpha,
                                           this->ld.gedm,
                                           this->ld.inl_index,
                                           fvnl_dalpha,
                                           cal_stress,
                                           svnl_dalpha);
    }
    else
    {
        const int nks = kv.nkstot;
        DeePKS_domain::cal_f_delta<std::complex<double>>(dm_k_new,
                                                         ucell,
                                                         ORB,
                                                         Test_Deepks::GridD,
                                                         ParaO,
                                                         nks,
                                                         kv.kvec_d,
                                                         this->ld.phialpha,
                                                         this->ld.gedm,
                                                         this->ld.inl_index,
                                                         fvnl_dalpha,
                                                         cal_stress,
                                                         svnl_dalpha);
    }
    DeePKS_domain::check_f_delta(ucell.nat, fvnl_dalpha, svnl_dalpha);

    this->compare_with_ref("F_delta.dat", "F_delta_ref.dat");
    this->compare_with_ref("stress_delta.dat", "stress_delta_ref.dat");
}

void test_deepks::compare_with_ref(const std::string f1, const std::string f2)
{
    this->total_check += 1;
    std::ifstream file1(f1.c_str());
    std::ifstream file2(f2.c_str());
    double test_thr = 1e-8;

    std::string word1;
    std::string word2;
    while (file1 >> word1)
    {
        file2 >> word2;
        if ((word1[0] - '0' >= 0 && word1[0] - '0' < 10) || word1[0] == '-')
        {
            double num1 = std::stof(word1);
            double num2 = std::stof(word2);
            if (std::abs(num1 - num2) > test_thr)
            {
                this->failed_check += 1;
                std::cout << "\e[1;31m [  FAILED  ] \e[0m" << f1.c_str() << " inconsistent!" << std::endl;
                return;
            }
        }
        else
        {
            if (word1 != word2)
            {
                this->failed_check += 1;
                return;
            }
        }
    }
    return;
}
