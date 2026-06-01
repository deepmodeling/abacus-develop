#include "write_H_terms.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_hs/write_HS.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/filename.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/operator_force_stress_utils.h"

#include <complex>
#include <tuple>

namespace ModuleIO
{

static void setup_veff_hcontainer(hamilt::HContainer<double>& hR,
                                  const UnitCell& ucell,
                                  const Grid_Driver& gd,
                                  const Parallel_Orbitals& pv,
                                  const std::vector<double>& orb_cutoff)
{
    const Parallel_Orbitals* paraV = hR.get_paraV();
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1 = 0, I1 = 0;
        ucell.iat2iait(iat1, &I1, &T1);

        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell.itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0 < orb_cutoff[T1] + orb_cutoff[T2])
            {
                hamilt::AtomPair<double> tmp(iat1, iat2, R_index, paraV);
                hR.insert_pair(tmp);
            }
        }
    }
    hR.allocate(nullptr, true);
}

static void gather_and_write(const std::string& prefix,
                             const std::string& label,
                             hamilt::HContainer<double>& hR,
                             const UnitCell& ucell,
                             const Parallel_Orbitals& pv,
                             const int nspin,
                             const int ispin,
                             const int istep,
                             const bool append,
                             const int* iat2iwt,
                             const int nat)
{
    const int nbasis = hR.get_nbasis();
#ifdef __MPI
    Parallel_Orbitals serialV;
    serialV.init(nbasis, nbasis, nbasis, pv.comm());
    serialV.set_serial(nbasis, nbasis);
    serialV.set_atomic_trace(iat2iwt, nat, nbasis);
    hamilt::HContainer<double> hr_serial(&serialV);
    hamilt::gatherParallels(hR, &hr_serial, 0);
    if (GlobalV::MY_RANK == 0)
#endif
    {
        std::string fname;
        if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
        {
            fname = PARAM.globalv.global_matrix_dir + hsr_gen_fname(prefix, ispin, append, istep);
        }
        else
        {
            fname = PARAM.globalv.global_out_dir + hsr_gen_fname(prefix, ispin, append, istep);
        }
#ifdef __MPI
        write_hcontainer_csr(fname, &ucell, 8, &hr_serial, istep, ispin, nspin, label);
#else
        write_hcontainer_csr(fname, &ucell, 8, &hR, istep, ispin, nspin, label);
#endif
    }
}

static void write_hk_common(hamilt::HContainer<double>& hR,
                            const std::string& prefix,
                            const UnitCell& ucell,
                            const Parallel_Orbitals& pv,
                            const K_Vectors& kv,
                            const int nspin,
                            const int istep,
                            const bool append,
                            const int* iat2iwt,
                            const int nat)
{
    const int nspin_k = (nspin == 2 ? 2 : 1);
    const int nks = kv.get_nks() / nspin_k;
    const int nlocal = PARAM.globalv.nlocal;
    const bool gamma_only = PARAM.globalv.gamma_only_local;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool out_app_flag = PARAM.inp.out_app_flag;

    for (int ik = 0; ik < nks; ++ik)
    {
        const ModuleBase::Vector3<double>& kvec_d = kv.kvec_d[ik];

        std::vector<std::complex<double>> hk_global(nlocal * nlocal, 0);
        hamilt::folding_HR(hR, hk_global.data(), kvec_d, nlocal, 0);

        const int out_label = 1;
        std::string fname = ModuleIO::filename_output(global_out_dir,
                                                      prefix,
                                                      "nao",
                                                      ik,
                                                      kv.ik2iktot,
                                                      nspin,
                                                      kv.get_nkstot(),
                                                      out_label,
                                                      out_app_flag,
                                                      gamma_only,
                                                      istep);
        std::cout << "hk term fname = " << fname << std::endl;
        ModuleIO::save_mat(istep,
                           hk_global.data(),
                           nlocal,
                           false,
                           8,
                           false,
                           out_app_flag,
                           fname,
                           pv,
                           GlobalV::DRANK);
    }
}

void write_h_t(const UnitCell& ucell,
               const Grid_Driver& gd,
               const Parallel_Orbitals& pv,
               const TwoCenterBundle& two_center_bundle,
               const LCAO_Orbitals& orb,
               const K_Vectors& kv,
               const int nspin,
               const int istep,
               const bool append,
               const int* iat2iwt,
               const int nat,
               const bool also_hk)
{
    ModuleBase::TITLE("ModuleIO", "write_h_t");
    ModuleBase::timer::start("ModuleIO", "write_h_t");

    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nspin_out = (nspin == 2 ? 2 : 1);

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));

        hamilt::EKinetic<hamilt::OperatorLCAO<double, double>>
            tmp_ekinetic(nullptr, kv.kvec_d, &hR_tmp, &ucell, orb_cutoff, &gd, two_center_bundle.kinetic_orb.get());
        tmp_ekinetic.contributeHR();

        gather_and_write("t", "T", hR_tmp, ucell, pv, nspin, ispin, istep, append, iat2iwt, nat);

        if (also_hk)
        {
            write_hk_common(hR_tmp, "tk", ucell, pv, kv, nspin, istep, append, iat2iwt, nat);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_t");
}

void write_h_vnl(const UnitCell& ucell,
                 const Grid_Driver& gd,
                 const Parallel_Orbitals& pv,
                 const TwoCenterBundle& two_center_bundle,
                 const LCAO_Orbitals& orb,
                 const K_Vectors& kv,
                 const int nspin,
                 const int istep,
                 const bool append,
                 const int* iat2iwt,
                 const int nat,
                 const bool also_hk)
{
    ModuleBase::TITLE("ModuleIO", "write_h_vnl");
    ModuleBase::timer::start("ModuleIO", "write_h_vnl");

    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nspin_out = (nspin == 2 ? 2 : 1);

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));

        hamilt::Nonlocal<hamilt::OperatorLCAO<double, double>> tmp_nonlocal(nullptr,
                                                                            kv.kvec_d,
                                                                            &hR_tmp,
                                                                            &ucell,
                                                                            orb_cutoff,
                                                                            &gd,
                                                                            two_center_bundle.overlap_orb_beta.get());
        tmp_nonlocal.contributeHR();

        gather_and_write("vnl", "V^NL", hR_tmp, ucell, pv, nspin, ispin, istep, append, iat2iwt, nat);

        if (also_hk)
        {
            write_hk_common(hR_tmp, "vnlk", ucell, pv, kv, nspin, istep, append, iat2iwt, nat);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_vnl");
}

void write_h_vl(const UnitCell& ucell,
                const Grid_Driver& gd,
                const Parallel_Orbitals& pv,
                const LCAO_Orbitals& orb,
                const elecstate::Potential* pot,
                const int nspin,
                const int istep,
                const bool append,
                const int* iat2iwt,
                const int nat,
    const K_Vectors& kv,
                const bool also_hk)
{
    ModuleBase::TITLE("ModuleIO", "write_h_vl");
    ModuleBase::timer::start("ModuleIO", "write_h_vl");

    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nspin_out = (nspin == 2 ? 2 : 1);

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);

        const double* v_local = pot->get_fixed_v(); // local pp, no Hxc
        ModuleGint::cal_gint_vl(v_local, &hR_tmp);

        gather_and_write("vl", "V^L", hR_tmp, ucell, pv, nspin, ispin, istep, append, iat2iwt, nat);

        if (also_hk)
        {
            write_hk_common(hR_tmp, "vlk", ucell, pv, kv, nspin, istep, append, iat2iwt, nat);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_vl");
}

void write_h_vh(const UnitCell& ucell,
                const Grid_Driver& gd,
                const Parallel_Orbitals& pv,
                const LCAO_Orbitals& orb,
                const Charge* chg,
                const ModulePW::PW_Basis* rho_basis,
                const int nspin,
                const int istep,
                const bool append,
                const int* iat2iwt,
                const int nat,
    const K_Vectors& kv,
                const bool also_hk)
{
    ModuleBase::TITLE("ModuleIO", "write_h_vh");
    ModuleBase::timer::start("ModuleIO", "write_h_vh");

    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nspin_out = (nspin == 2 ? 2 : 1);

    ModuleBase::matrix v_h
        = elecstate::H_Hartree_pw::v_hartree(ucell, const_cast<ModulePW::PW_Basis*>(rho_basis), nspin, chg->rho);

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);

        ModuleGint::cal_gint_vl(&v_h(ispin, 0), &hR_tmp);

        gather_and_write("vh", "V^H", hR_tmp, ucell, pv, nspin, ispin, istep, append, iat2iwt, nat);

        if (also_hk)
        {
            write_hk_common(hR_tmp, "vhk", ucell, pv, kv, nspin, istep, append, iat2iwt, nat);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_vh");
}

void write_h_vxc(const UnitCell& ucell,
                 const Grid_Driver& gd,
                 const Parallel_Orbitals& pv,
                 const LCAO_Orbitals& orb,
                 const Charge* chg,
                 const int nrxx,
                 const int nspin,
                 const int istep,
                 const bool append,
                 const int* iat2iwt,
                 const int nat,
    const K_Vectors& kv,
                 const bool also_hk)
{
    ModuleBase::TITLE("ModuleIO", "write_h_vxc");
    ModuleBase::timer::start("ModuleIO", "write_h_vxc");

    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nspin_out = (nspin == 2 ? 2 : 1);

    ModuleBase::matrix v_xc;
    double etxc, vtxc;
    std::tie(etxc, vtxc, v_xc) = XC_Functional::v_xc(nrxx, chg, &ucell);

    for (int ispin = 0; ispin < nspin_out; ispin++)
    {
        hamilt::HContainer<double> hR_tmp(const_cast<Parallel_Orbitals*>(&pv));
        setup_veff_hcontainer(hR_tmp, ucell, gd, pv, orb_cutoff);

        ModuleGint::cal_gint_vl(&v_xc(ispin, 0), &hR_tmp);

        gather_and_write("vxc", "V^XC", hR_tmp, ucell, pv, nspin, ispin, istep, append, iat2iwt, nat);

        if (also_hk)
        {
            write_hk_common(hR_tmp, "vxck", ucell, pv, kv, nspin, istep, append, iat2iwt, nat);
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_h_vxc");
}

} // namespace ModuleIO
