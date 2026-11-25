#include "write_dmr.h"
#include "source_io/module_hs/hs_sparse_io_detail.h"
#include "source_io/module_output/spin_tag.h"
#include "source_io/module_output/write_notice.h"

#include "source_base/global_variable.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"

namespace ModuleIO
{
std::string dmr_gen_fname(const int out_type, const int ispin, const bool numbered, const int istep)
{
    const std::string spin_suffix = "s" + std::to_string(ispin + 1);
    if (out_type == 1)
    {
        // spa stands for sparse
        return detail::make_sparse_R_filename("dmr", spin_suffix, numbered, istep, ".csr");
    }
    else if (out_type == 2)
    {
        return detail::make_sparse_R_filename("dmr", spin_suffix, false, istep, ".npz");
    }
    ModuleBase::WARNING("write_dmr", "the output type of density matrix DM(R) should be csr or npz.");
    return "dmr.csr";
}

void write_dmr_csr(std::string& fname,
                   const UnitCell *ucell,
                   const int precision,
                   hamilt::HContainer<double>* dm_serial,
                   const int istep,
                   const int ispin,
                   const int nspin,
                   const bool append)
{
    // write the head: ION step number, basis number and R loop number

    std::ofstream ofs;

    // The first written step of a run (or a single-point SCF) truncates the
    // fixed file so a rerun does not keep stale records; later steps append
    // only when the caller asks for it (mode 1).
    const bool append_now = append && (istep > 0);
    if (append_now)
    {
        ofs.open(fname, std::ios::app);
    }
    else
    {
        ofs.open(fname);
    }
    ModuleIO::detail::check_output_file_open(ofs, fname, "ModuleIO::write_dmr_csr");

    const std::string description = " # print density matrix in real space DM(R)";
    const ModuleIO::detail::CsrTextHeader header = {
        istep + 1,
        description,
        nspin,
        ispin + 1,
        dm_serial->get_nbasis(),
        static_cast<int>(dm_serial->size_R_loop()),
        ucell,
        std::string(),
        false,
        0.0,
        false};
    ModuleIO::detail::write_csr_text_header(ofs, header);

    // write HR_serial to ofs
    const double sparse_threshold = 1e-10;
    hamilt::Output_HContainer<double> dmr(dm_serial, ofs, sparse_threshold, precision);
    dmr.write();
    ofs.close();
}

void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr,
               const UnitCell* ucell,
               const int precision,
               const Parallel_Blkcyclic& paraV,
               const int out_app_flag,
               const int* iat2iwt,
               const int nat,
               const std::string& dmr_dir,
               const int istep,
               std::ofstream& ofs_running)
{
    const int nspin = dmr.size();
    assert(nspin > 0);
    // Mode 1 appends one frame per ionic step; the first step (istep == 0)
    // truncates the fixed file so a rerun does not keep stale records. Only
    // mode 2 embeds the step in the name.
    const bool append = (out_app_flag == 1) && (istep > 0);
    const bool numbered = (out_app_flag == 2);
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        const int nbasis = dmr[ispin]->get_nbasis();

        // gather the parallel matrix to serial matrix
#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.para_2d().set_serial(nbasis, nbasis);
        serialV.set_atomic_trace(iat2iwt, nat, nbasis);
        hamilt::HContainer<double> dm_serial(&serialV);
        hamilt::gatherParallels(*dmr[ispin], &dm_serial, 0);
#else
        hamilt::HContainer<double> dm_serial(*dmr[ispin]);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            // out_type = 1, csr format;
            // out_type = 2, npz format (currently not support)
            const int out_type = 1;
            std::string fname = dmr_dir + dmr_gen_fname(out_type, ispin, numbered, istep);
            write_dmr_csr(fname, ucell, precision, &dm_serial, istep, ispin, nspin, append);
            ModuleIO::write_matrix_notice(ofs_running, "DM(R)", ispin, nspin, fname);
        }
    }
}

} // namespace ModuleIO
