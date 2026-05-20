#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "write_dH.h"

namespace ModuleIO
{

bool write_dH_vh(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh");

    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dVh_x(&pv);
        hamilt::HContainer<double> dVh_y(&pv);
        hamilt::HContainer<double> dVh_z(&pv);
        dVh_x.set_zero();
        dVh_y.set_zero();
        dVh_z.set_zero();

        const int nbasis = dVh_x.get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dVh_x_s(&serialV);
        hamilt::HContainer<double> dVh_y_s(&serialV);
        hamilt::HContainer<double> dVh_z_s(&serialV);
        hamilt::gatherParallels(dVh_x, &dVh_x_s, 0);
        hamilt::gatherParallels(dVh_y, &dVh_y_s, 0);
        hamilt::gatherParallels(dVh_z, &dVh_z_s, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dvhrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dvhry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dvhrz", ispin, params.append, params.istep);

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
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVh_x_s, params.istep, ispin, nspin, "dV^H");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVh_y_s, params.istep, ispin, nspin, "dV^H");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVh_z_s, params.istep, ispin, nspin, "dV^H");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVh_x, params.istep, ispin, nspin, "dV^H");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVh_y, params.istep, ispin, nspin, "dV^H");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVh_z, params.istep, ispin, nspin, "dV^H");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vh");
    return true;
}

} // namespace ModuleIO
