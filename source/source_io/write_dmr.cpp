#include "write_dmr.h"

#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_pw/module_pwdft/global.h"

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

void write_dmr_csr(std::string& fname, hamilt::HContainer<double>* dm_serial, const int istep)
{
    // write the head: ION step number, basis number and R loop number

	std::ofstream ofs;

    // mohan update 2025-05-26
	if(istep<=0)
	{
		ofs.open(fname);
	}
	else if(istep>0)
	{
		ofs.open(fname, std::ios::app);
	}

    ofs << "IONIC_STEP: " << istep+1 << std::endl;
    ofs << "Matrix Dimension of DM(R): " << dm_serial->get_nbasis() << std::endl;
    ofs << "Matrix number of DM(R): " << dm_serial->size_R_loop() << std::endl;

    // write HR_serial to ofs
    const double sparse_threshold = 1e-10;
    const int precision = 8;
    hamilt::Output_HContainer<double> out_dmr(dm_serial, ofs, sparse_threshold, precision);
    out_dmr.write();
    ofs.close();
}

void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr,
               const Parallel_2D& paraV,
               const bool append,
               const int* iat2iwt,
               const int nat,
               const int istep)
{
    for (int ispin = 0; ispin < dmr.size(); ispin++)
	{
		const int nbasis = dmr[ispin]->get_nbasis();
		// gather the parallel matrix to serial matrix
#ifdef __MPI
		Parallel_Orbitals serialV;
		serialV.init(nbasis, nbasis, nbasis, paraV.comm());
		serialV.set_serial(nbasis, nbasis);
		serialV.set_atomic_trace(iat2iwt, nat, nbasis);
		hamilt::HContainer<double> dm_serial(&serialV);
		hamilt::gatherParallels(*dmr[ispin], &dm_serial, 0);
#else
		hamilt::HContainer<double> dm_serial(*dmr[ispin]);
#endif
		if (GlobalV::MY_RANK == 0)
		{
			std::string fname = PARAM.globalv.global_out_dir + dmr_gen_fname(1, ispin, append, istep);
			write_dmr_csr(fname, &dm_serial, istep);
		}
	}
}

} // namespace ModuleIO
