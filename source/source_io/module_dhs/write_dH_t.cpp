#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/operator_force_stress_utils.h"
#include "write_dH.h"

namespace ModuleIO
{

bool write_dH_t(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_t");
    ModuleBase::timer::start("ModuleIO", "write_dH_t");

    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const LCAO_Orbitals& orb = *params.orb;

    const int nat = ucell.nat;
    const int nspin = params.nspin;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dT_x(&pv);
        hamilt::HContainer<double> dT_y(&pv);
        hamilt::HContainer<double> dT_z(&pv);

        hamilt::EKinetic<hamilt::OperatorLCAO<double, double>> tmp_ekinetic(nullptr,
                                                                            params.kv->kvec_d,
                                                                            nullptr,
                                                                            &ucell,
                                                                            orb_cutoff,
                                                                            &gd,
                                                                            two_center_bundle.kinetic_orb.get());

        tmp_ekinetic.cal_dH(&dT_x, &dT_y, &dT_z);

        const int nbasis = dT_x.get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dT_x_serial(&serialV);
        hamilt::HContainer<double> dT_y_serial(&serialV);
        hamilt::HContainer<double> dT_z_serial(&serialV);
        hamilt::gatherParallels(dT_x, &dT_x_serial, 0);
        hamilt::gatherParallels(dT_y, &dT_y_serial, 0);
        hamilt::gatherParallels(dT_z, &dT_z_serial, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dtrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dtry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dtrz", ispin, params.append, params.istep);

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
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dT_x_serial, params.istep, ispin, nspin, "dT");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dT_y_serial, params.istep, ispin, nspin, "dT");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dT_z_serial, params.istep, ispin, nspin, "dT");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dT_x, params.istep, ispin, nspin, "dT");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dT_y, params.istep, ispin, nspin, "dT");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dT_z, params.istep, ispin, nspin, "dT");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_t");
    return true;
}

} // namespace ModuleIO
