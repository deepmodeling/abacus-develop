#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "write_dH.h"

namespace ModuleIO
{

bool write_dH_vxc(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vxc");
    ModuleBase::timer::start("ModuleIO", "write_dH_vxc");

    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dVxc_x(&pv);
        hamilt::HContainer<double> dVxc_y(&pv);
        hamilt::HContainer<double> dVxc_z(&pv);
        dVxc_x.set_zero();
        dVxc_y.set_zero();
        dVxc_z.set_zero();

        const int nbasis = dVxc_x.get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dVxc_x_s(&serialV);
        hamilt::HContainer<double> dVxc_y_s(&serialV);
        hamilt::HContainer<double> dVxc_z_s(&serialV);
        hamilt::gatherParallels(dVxc_x, &dVxc_x_s, 0);
        hamilt::gatherParallels(dVxc_y, &dVxc_y_s, 0);
        hamilt::gatherParallels(dVxc_z, &dVxc_z_s, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dvxcrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dvxcry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dvxcrz", ispin, params.append, params.istep);

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
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVxc_x_s, params.istep, ispin, nspin, "dV^XC");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVxc_y_s, params.istep, ispin, nspin, "dV^XC");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVxc_z_s, params.istep, ispin, nspin, "dV^XC");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVxc_x, params.istep, ispin, nspin, "dV^XC");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVxc_y, params.istep, ispin, nspin, "dV^XC");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVxc_z, params.istep, ispin, nspin, "dV^XC");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vxc");
    return true;
}

} // namespace ModuleIO
