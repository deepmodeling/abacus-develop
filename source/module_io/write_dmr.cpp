#include <iostream>

#include "write_dmr.h"
#include "module_hamilt_lcao/module_hcontainer/output_hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleIO
{
std::string dmr_gen_fname(const int out_type, const bool sparse, const int ispin, const bool append, const int istep)
{
    std::string fname="dmr.csr";
    if (out_type == 1)
    {
        if (sparse)
        {
            if (!append && istep >=0)
            {
                fname = std::to_string(istep + 1) + "_data-DMR-sparse_SPIN" + std::to_string(ispin) + ".csr";
            }
            else
            {
                fname = "data-DMR-sparse_SPIN" + std::to_string(ispin) + ".csr";
            }
        }
        else
        {
            fname = "SPIN" + std::to_string(ispin + 1) + "_DM";
        }
    }
    else if (out_type == 2 && ! sparse)
    {
        if (sparse)
        {
            fname = "output_DM" + std::to_string(ispin) + ".npz";
        }
        else
        {
            ModuleBase::WARNING("write_dmr","do not support the output of dense DMR in npz type.");
        }
    }
    else
    {
        ModuleBase::WARNING("write_dmr","the output type of DMR should be npz or csr.");
    }
    return fname;
}


void write_dmr_csr(std::string& fname, const hamilt::HContainer<double>& dm, const int istep)
{
    // gather dm to dm_serial
    Parallel_Orbitals serialV;
    serialV.set_serial(GlobalV::NLOCAL, GlobalV::NLOCAL);
    serialV.set_atomic_trace(GlobalC::ucell.get_iat2iwt(), GlobalC::ucell.nat, GlobalV::NLOCAL);
    hamilt::HContainer<double>* dm_serial; 
#ifdef __MPI    
    if (GlobalV::MY_RANK == 0)
    {
        dm_serial = new hamilt::HContainer<double>(&serialV);
    }
    hamilt::gatherParallels(dm, dm_serial, 0);
#else
    dm_serial = new hamilt::HContainer<double>(dm);
#endif

    // write the head: ION step number, basis number and R loop number
    std::ofstream ofs(fname, std::ios::app);
    ofs << "STEP: " << istep << std::endl;
    ofs << "Matrix Dimension of DM(R): " << GlobalV::NLOCAL << std::endl;
    ofs << "Matrix number of DM(R): " << dm_serial->size_R_loop() << std::endl;

    // write HR_serial to ofs
    double sparse_threshold = 1e-10;
    int precision = 8;
    hamilt::Output_HContainer<double> out_dm(dm_serial, &serialV, ofs, sparse_threshold, precision);
    out_dm.write();
    ofs.close();
    delete dm_serial;
}

void write_dmr(
    const hamilt::HContainer<double>& dm,
    const int out_type,
    const bool sparse,
    const int ispin,
    const bool append,
    const int istep, 
    const Parallel_Orbitals& pv
    )
{
    if (out_type != 1 && out_type != 2)
    {
        ModuleBase::WARNING("write_dmr","Only support output npz or csr type now.");
        return;
    }

    std::string fname = GlobalV::global_out_dir + dmr_gen_fname(out_type, sparse, ispin, append, istep);

    if (out_type == 1 && sparse) // for out_dm1
    {
        write_dmr_csr(fname, dm, istep);
    }
    else if (out_type == 2 && sparse) // for out_dm_npz
    {
        //output_mat_npz(fname,dm);
    }
    else
    {
        ModuleBase::WARNING("write_dmr","Only support to output the sparse DMR in npz or csr type now.");
        return;
    }
}

}