#include "write_dmr.h"
#include "source_io/module_output/spin_tag.h"

#include "source_base/global_variable.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"
#include "source_cell/ucell_io.h"

namespace ModuleIO
{
std::string dmr_gen_fname(const int out_type, const int ispin, const bool numbered, const int istep)
{
    std::string fname = "dmr.csr";
    if (out_type == 1)
    {
        if (numbered && istep >= 0)
        {
            // spa stands for sparse
            fname = "dmrs" + std::to_string(ispin+1) + "g" + std::to_string(istep + 1) + "_nao.csr";
        }
        else
        {
            fname = "dmrs" + std::to_string(ispin+1) + "_nao.csr";
        }
    }
    else if (out_type == 2)
    {
        fname = "dmrs" + std::to_string(ispin+1) + "_nao.npz";
    }
    else
    {
        ModuleBase::WARNING("write_dmr", "the output type of density matrix DM(R) should be csr or npz.");
    }
    return fname;
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


    ofs << "Ionic_Step " << istep+1 << std::endl;
    ofs << " # print density matrix in real space DM(R)" << std::endl;
    ofs << " " << nspin << " # number of spin directions" << std::endl;
    ofs << " " << ispin+1 << " # spin index" << std::endl;
    ofs << " " << dm_serial->get_nbasis()
        << " # number of localized basis" << std::endl;
    ofs << " " << dm_serial->size_R_loop()
        << " # number of Bravais lattice vector R" << std::endl;
    ofs << std::endl;

    // write ucell
    ModuleIO::UcellIO::write_ucell(ofs, ucell);
    ofs << std::endl;

    // write HR_serial to ofs
    const double sparse_threshold = 1e-10;
    hamilt::Output_HContainer<double> dmr(dm_serial, ofs, sparse_threshold, precision);
    dmr.write();
    ofs.close();
}

void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr,
               const UnitCell* ucell,
               const int precision,
               const BlockCyclicDistribution& paraV,
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
            const std::string spin_tag = ModuleIO::make_spin_tag(ispin, nspin);
            ofs_running << " Write DM(R)" << spin_tag << " matrix in NAO basis to file: " << fname << std::endl;
        }
    }
}

} // namespace ModuleIO
