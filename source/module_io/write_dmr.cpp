#include "write_dmr.h"

#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hamilt_lcao/module_hcontainer/output_hcontainer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include <iostream>

namespace ModuleIO
{
std::string dmr_gen_fname(const int out_type, const int ispin, const bool append, const int istep)
{
    std::string fname = "dmr.csr";
    if (out_type == 1)
    {
        if (!append && istep >= 0)
        {
            fname = std::to_string(istep + 1) + "_data-DMR-sparse_SPIN" + std::to_string(ispin) + ".csr";
        }
        else
        {
            fname = "data-DMR-sparse_SPIN" + std::to_string(ispin) + ".csr";
        }
    }
    else if (out_type == 2)
    {
        fname = "output_DM" + std::to_string(ispin) + ".npz";
    }
    else
    {
        ModuleBase::WARNING("write_dmr", "the output type of DMR should be npz or csr.");
    }
    return fname;
}

void write_dmr_csr(std::string& fname, hamilt::HContainer<double>* dm_serial, const int istep)
{
    // write the head: ION step number, basis number and R loop number
    std::ofstream ofs(fname, std::ios::app);
    ofs << "STEP: " << istep << std::endl;
    ofs << "Matrix Dimension of DM(R): " << GlobalV::NLOCAL << std::endl;
    ofs << "Matrix number of DM(R): " << dm_serial->size_R_loop() << std::endl;

    // write HR_serial to ofs
    double sparse_threshold = 1e-10;
    int precision = 8;
    hamilt::Output_HContainer<double> out_dm(dm_serial, GlobalV::NLOCAL, GlobalV::NLOCAL, ofs, sparse_threshold, precision);
    out_dm.write();
    ofs.close();
}

void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr, const bool out_csr, const bool out_npz,
               const bool append, const int istep)
{
    if (!out_csr && !out_npz)
    {
        ModuleBase::WARNING("write_dmr", "the output type of DMR should be npz or csr.");
        return;
    }

    for (int ispin = 0; ispin < dmr.size(); ispin++)
    {
        // gather the parallel matrix to serial matrix
        hamilt::HContainer<double>* dm_serial;
#ifdef __MPI
        Parallel_Orbitals serialV;
        serialV.set_serial(GlobalV::NLOCAL, GlobalV::NLOCAL);
        serialV.set_atomic_trace(GlobalC::ucell.get_iat2iwt(), GlobalC::ucell.nat, GlobalV::NLOCAL);
        dm_serial = new hamilt::HContainer<double>(&serialV);
        hamilt::gatherParallels(*dmr[ispin], dm_serial, 0);
#else
        dm_serial = new hamilt::HContainer<double>(*dmr[ispin]);
#endif

        if(out_csr)
        {
            if (GlobalV::MY_RANK == 0)
            {
                std::string fname = GlobalV::global_out_dir + dmr_gen_fname(1, ispin, append, istep);
                write_dmr_csr(fname, dm_serial, istep);
            }
        }

        if(out_npz)
        {
            // std::string fname = GlobalV::global_out_dir + dmr_gen_fname(2, ispin, append, istep);
            // write_dmr_npz(fname, *dm_serial, istep);
        }
        delete dm_serial;
    }
}

} // namespace ModuleIO