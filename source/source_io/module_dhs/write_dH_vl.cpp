#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_lcao/module_gint/phi_operator.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/operator_force_stress_utils.h"
#include "write_dH.h"

#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleIO
{

static void write_dH_vl_pulay(const WriteDHParams& params,
                              int ispin,
                              hamilt::HContainer<double>& dVl_x,
                              hamilt::HContainer<double>& dVl_y,
                              hamilt::HContainer<double>& dVl_z)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const LCAO_Orbitals& orb = *params.orb;
    const std::vector<double>& orb_cutoff = orb.cutoffs();
    const int nat = ucell.nat;
    const int npol = ucell.get_npol();

#pragma omp parallel for schedule(dynamic)
    for (int iat = 0; iat < nat; iat++)
    {
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        const ModuleBase::Vector3<double>& tau1 = ucell.atoms[T1].tau[I1];
        const Atom* atom1 = &ucell.atoms[T1];

        AdjacentAtomInfo adjs;
        gd.Find_atom(ucell, tau1, T1, I1, &adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ad++)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell.itia2iat(T2, I2);
            const ModuleBase::Vector3<double>& tau2 = adjs.adjacent_tau[ad];
            const Atom* atom2 = &ucell.atoms[T2];

            ModuleBase::Vector3<double> dtau = tau2 - tau1;
            double dist = dtau.norm() * ucell.lat0;
            double cutoff = orb_cutoff[T1] + orb_cutoff[T2];

            if (dist >= cutoff)
                continue;

            ModuleBase::Vector3<int> R_idx(adjs.box[ad].x, adjs.box[ad].y, adjs.box[ad].z);

            hamilt::BaseMatrix<double>* mtx_x = dVl_x.find_matrix(iat, iat2, R_idx.x, R_idx.y, R_idx.z);
            hamilt::BaseMatrix<double>* mtx_y = dVl_y.find_matrix(iat, iat2, R_idx.x, R_idx.y, R_idx.z);
            hamilt::BaseMatrix<double>* mtx_z = dVl_z.find_matrix(iat, iat2, R_idx.x, R_idx.y, R_idx.z);

            if (!mtx_x || !mtx_y || !mtx_z)
                continue;

            double* ptr_x = mtx_x->get_pointer();
            double* ptr_y = mtx_y->get_pointer();
            double* ptr_z = mtx_z->get_pointer();
            const int row_sz = mtx_x->get_row_size();
            const int col_sz = mtx_x->get_col_size();

            for (int ii = 0; ii < atom1->nw * npol; ii++)
            {
                const int iw1_all = ucell.itiaiw2iwt(T1, I1, ii / npol);
                if (pv.global2local_row(iw1_all) < 0)
                    continue;
                const int iw1_row = pv.global2local_row(iw1_all);

                for (int jj = 0; jj < atom2->nw * npol; jj++)
                {
                    const int iw2_all = ucell.itiaiw2iwt(T2, I2, jj / npol);
                    if (pv.global2local_col(iw2_all) < 0)
                        continue;
                    const int iw2_col = pv.global2local_col(iw2_all);

                    int idx = iw1_row * col_sz + iw2_col;

                    // Contribution from atom1's orbitals: -<grad phi_iw1|V^L|phi_iw2>
                    // This is the Pulay term which needs grid integration
                    // Placeholder: store zero for now
                    ptr_x[idx] += 0.0;
                    ptr_y[idx] += 0.0;
                    ptr_z[idx] += 0.0;
                }
            }
        }
    }
}

bool write_dH_vl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vl");

    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> dVl_x(&pv);
        hamilt::HContainer<double> dVl_y(&pv);
        hamilt::HContainer<double> dVl_z(&pv);
        dVl_x.set_zero();
        dVl_y.set_zero();
        dVl_z.set_zero();

        // Compute Pulay term via grid integration
        write_dH_vl_pulay(params, ispin, dVl_x, dVl_y, dVl_z);

        const int nbasis = dVl_x.get_nbasis();

        // Write to CSR files
#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.init(nbasis, nbasis, nbasis, pv.comm());
        serialV.set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(params.iat2iwt, nat, nbasis);

        hamilt::HContainer<double> dVl_x_s(&serialV);
        hamilt::HContainer<double> dVl_y_s(&serialV);
        hamilt::HContainer<double> dVl_z_s(&serialV);
        hamilt::gatherParallels(dVl_x, &dVl_x_s, 0);
        hamilt::gatherParallels(dVl_y, &dVl_y_s, 0);
        hamilt::gatherParallels(dVl_z, &dVl_z_s, 0);

        if (GlobalV::MY_RANK == 0)
#endif
        {
            std::string fname_x = ModuleIO::dhr_gen_fname("dvlrx", ispin, params.append, params.istep);
            std::string fname_y = ModuleIO::dhr_gen_fname("dvlry", ispin, params.append, params.istep);
            std::string fname_z = ModuleIO::dhr_gen_fname("dvlrz", ispin, params.append, params.istep);

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
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVl_x_s, params.istep, ispin, nspin, "dV^L");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVl_y_s, params.istep, ispin, nspin, "dV^L");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVl_z_s, params.istep, ispin, nspin, "dV^L");
#else
            ModuleIO::write_hcontainer_csr(fname_x, &ucell, 8, &dVl_x, params.istep, ispin, nspin, "dV^L");
            ModuleIO::write_hcontainer_csr(fname_y, &ucell, 8, &dVl_y, params.istep, ispin, nspin, "dV^L");
            ModuleIO::write_hcontainer_csr(fname_z, &ucell, 8, &dVl_z, params.istep, ispin, nspin, "dV^L");
#endif
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vl");
    return true;
}

} // namespace ModuleIO
