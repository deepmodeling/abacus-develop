#include "write_dH.h"

#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_io/module_hs/write_HS.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"

#include <algorithm>
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
    std::array<std::vector<hamilt::HContainer<double>*>, 3>& g,
    const std::vector<int>& atom_filter)
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

    const bool filter_atoms = !atom_filter.empty();
    if (filter_atoms)
        for (int idx : atom_filter)
            if (idx < 0 || idx >= nat)
                ModuleBase::WARNING("write_dh_perI",
                    "atom index " + std::to_string(idx + 1) + " (1-based) is out of range [1, "
                    + std::to_string(nat) + "] and will be skipped");
    for (int iat = 0; iat < nat; ++iat)
    {
        if (filter_atoms && std::find(atom_filter.begin(), atom_filter.end(), iat) == atom_filter.end())
            continue;
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

void write_dH_components(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_components");
    ModuleBase::timer::start("ModuleIO", "write_dH_components");

    // nspin=4 (noncollinear) is not supported: needs complex spinor blocks (HContainer<std::complex<double>>)
    // plus noncollinear Gint kernels that do not exist for the dvlocal/drho paths. 
    if (PARAM.inp.nspin == 4)
    {
        ModuleBase::WARNING_QUIT("write_dH_components",
                                 "dH/dR component output (out_mat_dh_*) is not supported for "
                                 "nspin=4 (noncollinear) yet; only nspin=1 and nspin=2.");
    }

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out dH/dR components#                       |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    if (PARAM.inp.out_mat_dh[0])
    {
        write_dH_sum(params);
    }

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

#ifdef __EXX
    if (PARAM.inp.out_mat_dh_exx[0])
    {
        write_dH_exx(params);
    }
#endif

    ModuleBase::timer::end("ModuleIO", "write_dH_components");
}

} // namespace ModuleIO
