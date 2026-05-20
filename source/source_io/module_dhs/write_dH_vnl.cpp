#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/operator_force_stress_utils.h"
#include "write_dH.h"

namespace ModuleIO
{

bool write_dH_vnl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vnl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vnl");

    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const int nat = ucell.nat;
    const int nspin = params.nspin;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dVnl_x(&pv);
        hamilt::HContainer<double> dVnl_y(&pv);
        hamilt::HContainer<double> dVnl_z(&pv);

        hamilt::Nonlocal<hamilt::OperatorLCAO<double, double>> tmp_nonlocal(nullptr,
                                                                            params.kv->kvec_d,
                                                                            nullptr,
                                                                            &ucell,
                                                                            orb_cutoff,
                                                                            &gd,
                                                                            two_center_bundle.overlap_orb_beta.get());

        tmp_nonlocal.cal_dH(&dVnl_x, &dVnl_y, &dVnl_z);

        const int nbasis = dVnl_x.get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dVnl_x_s(&serialV);
        hamilt::HContainer<double> dVnl_y_s(&serialV);
        hamilt::HContainer<double> dVnl_z_s(&serialV);
        hamilt::gatherParallels(dVnl_x, &dVnl_x_s, 0);
        hamilt::gatherParallels(dVnl_y, &dVnl_y_s, 0);
        hamilt::gatherParallels(dVnl_z, &dVnl_z_s, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dvnlrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dvnlry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dvnlrz", ispin, params.append, params.istep);

            if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
            {
                fname_x = PARAM.globalv.global_matrix_dir + fname_x;
                fname_y = PARAM.globalv.global_matrix_dir + fname_y;
                fname_z = PARAM.globalv.global_matrix_dir + fname_z;
            }
            else
            {
                fname_x = PARAM.globalv.global_out_dir + fname_x;
                fname_y = PARAM.globalv.global_out_dir + fname_y;
                fname_z = PARAM.globalv.global_out_dir + fname_z;
            }

#ifdef __MPI
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVnl_x_s, params.istep, ispin, nspin, "dV^NL");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVnl_y_s, params.istep, ispin, nspin, "dV^NL");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVnl_z_s, params.istep, ispin, nspin, "dV^NL");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVnl_x, params.istep, ispin, nspin, "dV^NL");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVnl_y, params.istep, ispin, nspin, "dV^NL");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVnl_z, params.istep, ispin, nspin, "dV^NL");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vnl");
    return true;
}

} // namespace ModuleIO
