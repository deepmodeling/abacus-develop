#include "write_dH.h"

#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_io/module_hs/write_HS.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"

#include <complex>
#include <fstream>
#include <iomanip>
#include <string>

namespace ModuleIO
{

void write_dh_perI(WriteDHParams& params,
                   int ispin,
                   const std::string& rprefix,
                   const std::string& kprefix,
                   const std::string& label,
    std::array<std::vector<hamilt::HContainer<double>*>, 3>& g)
{
    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = params.nat;
    const int nspin = params.nspin;
    const int nbasis = g[0][0]->get_nbasis();

    const char dirc[3] = { 'x', 'y', 'z' };

    // k-space (dense, folded like H(k)) parameters
    const int nspin_k = (nspin == 2 ? 2 : 1);
    const int nks = params.kv->get_nks() / nspin_k;
    const int nlocal = PARAM.globalv.nlocal;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool out_app_flag = PARAM.inp.out_app_flag;
    const std::string r_dir
        = (PARAM.inp.calculation == "md" && !out_app_flag) ? PARAM.globalv.global_matrix_dir : global_out_dir;

#ifdef __MPI
    Parallel_Orbitals serialV;
    serialV.init(nbasis, nbasis, nbasis, pv.comm());
    serialV.set_serial(nbasis, nbasis);
    serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);
#endif

    for (int iat = 0; iat < nat; ++iat)
    {
        for (int d = 0; d < 3; ++d)
        {
            hamilt::HContainer<double>* hR = g[d][iat];
            const std::string tag = std::string(1, dirc[d]) + "_iat" + std::to_string(iat + 1);

            // ---- real space dH(R), CSR (only when also_dhR; dH(k) below is always written) ----
            if (params.also_dhR)
            {
#ifdef __MPI
            hamilt::HContainer<double> hR_s(&serialV);
            hamilt::gatherParallels(*hR, &hR_s, 0);
            if (GlobalV::MY_RANK == 0)
#endif
            {
                std::string fr = r_dir + ModuleIO::dhr_gen_fname(rprefix + tag, ispin, params.append, params.istep);
#ifdef __MPI
                ModuleIO::write_hcontainer_csr(fr, &ucell, 8, &hR_s, params.istep, ispin, nspin, label);
#else
                ModuleIO::write_hcontainer_csr(fr, &ucell, 8, hR, params.istep, ispin, nspin, label);
#endif
            }
            }

            // ---- k space dH(k), dense (folded like H(k), comparable to *_nao.txt) ----
            // build the filename directly (filename_output only accepts a fixed property set)
#ifdef __MPI
            const bool col_major = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver);
            const size_t hk_size = static_cast<size_t>(pv.get_row_size()) * pv.get_col_size();
#else
            const size_t hk_size = static_cast<size_t>(nlocal) * nlocal;
#endif
            for (int ik = 0; ik < nks; ++ik)
            {
                std::vector<std::complex<double>> hk(hk_size, 0);
#ifdef __MPI
                if (col_major)
                    hamilt::folding_HR(*hR, hk.data(), params.kv->kvec_d[ik], pv.get_row_size(), 1);
                else
                    hamilt::folding_HR(*hR, hk.data(), params.kv->kvec_d[ik], pv.get_col_size(), 0);
#else
                hamilt::folding_HR(*hR, hk.data(), params.kv->kvec_d[ik], nlocal, 0);
#endif
                std::string fk = global_out_dir + kprefix + tag;
                if (nks > 1)
                {
                    fk += "_ik" + std::to_string(params.kv->ik2iktot[ik]);
                }
                fk += "_nao.txt";
                ModuleIO::save_mat(params.istep,
                                   hk.data(),
                                   nlocal,
                                   false,
                                   8,
                                   false,
                                   out_app_flag,
                                   fk,
                                   pv,
                                   GlobalV::DRANK);
            }
        }
    }
}

bool any_dh_term_enabled()
{
    return PARAM.inp.out_mat_dh_t[0] || PARAM.inp.out_mat_dh_vl[0] || PARAM.inp.out_mat_dh_vnl[0]
           || PARAM.inp.out_mat_dh_vh[0] || PARAM.inp.out_mat_dh_vxc[0];
}

void write_dH_components(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_components");
    ModuleBase::timer::start("ModuleIO", "write_dH_components");

    if (!any_dh_term_enabled())
    {
        ModuleBase::timer::end("ModuleIO", "write_dH_components");
        return;
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out dH/dR components#                       |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    if (PARAM.inp.out_mat_dh_t[0])
    {
        write_dH_t(params);
    }

    if (PARAM.inp.out_mat_dh_vnl[0])
    {
        write_dH_vnl(params);
    }

    if (PARAM.inp.out_mat_dh_vl[0])
    {
        write_dH_vl(params);
    }

    if (PARAM.inp.out_mat_dh_vh[0])
    {
        write_dH_vh(params);
        write_dH_vh_pulay(params);
    }

    if (PARAM.inp.out_mat_dh_vxc[0])
    {
        write_dH_vxc(params);
        write_dH_vxc_pulay(params);
    }

    if (PARAM.inp.out_mat_dh[0] && any_dh_term_enabled())
    {
        write_dH_sum(params);
    }

#ifdef __EXX
    if (PARAM.inp.out_mat_dh_exx[0])
    {
        write_dH_exx(params);
    }
#endif

    ModuleBase::timer::end("ModuleIO", "write_dH_components");
}

bool write_dH_sum(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_sum");
    ModuleBase::timer::start("ModuleIO", "write_dH_sum");

    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dH_sum_x(&pv);
        hamilt::HContainer<double> dH_sum_y(&pv);
        hamilt::HContainer<double> dH_sum_z(&pv);
        dH_sum_x.set_zero();
        dH_sum_y.set_zero();
        dH_sum_z.set_zero();

        if (dH_sum_x.size_atom_pairs() == 0)
        {
            continue;
        }

        const int nbasis = dH_sum_x.get_nbasis();

#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dH_x_s(&serialV);
        hamilt::HContainer<double> dH_y_s(&serialV);
        hamilt::HContainer<double> dH_z_s(&serialV);
        hamilt::gatherParallels(dH_sum_x, &dH_x_s, 0);
        hamilt::gatherParallels(dH_sum_y, &dH_y_s, 0);
        hamilt::gatherParallels(dH_sum_z, &dH_z_s, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dhrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dhry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dhrz", ispin, params.append, params.istep);

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
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dH_x_s, params.istep, ispin, nspin, "dH_sum");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dH_y_s, params.istep, ispin, nspin, "dH_sum");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dH_z_s, params.istep, ispin, nspin, "dH_sum");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dH_sum_x, params.istep, ispin, nspin, "dH_sum");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dH_sum_y, params.istep, ispin, nspin, "dH_sum");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dH_sum_z, params.istep, ispin, nspin, "dH_sum");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_sum");
    return true;
}

} // namespace ModuleIO
