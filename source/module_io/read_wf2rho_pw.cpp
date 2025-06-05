#include "read_wf2rho_pw.h"

#include "read_wfc_pw.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/kernels/elecstate_op.h"
#include "module_io/filename.h"

void ModuleIO::read_wf2rho_pw(
		const ModulePW::PW_Basis_K* pw_wfc,
		ModuleSymmetry::Symmetry& symm,
		Charge& chg,
        const std::string &readin_dir,
		const int kpar,
		const int my_pool,
		const int my_rank,
        const int nproc_in_pool,
        const int rank_in_pool,
		const int nbands,
		const int nspin,
		const int npol,
		const int nkstot,
		const std::vector<int> &ik2iktot,
		const std::vector<int> &isk,
		std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "read_wf2rho_pw");
    ModuleBase::timer::tick("ModuleIO", "read_wf2rho_pw");

	ofs_running << " READING WAVE FUNCTIONS" << std::endl;
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_running << " |                                            "
		"                        |" << std::endl;
	ofs_running << " | Reading electronic wave functions in plane wave basis set and      |" << std::endl;
	ofs_running << " | evaluate charge density based on these wave functions              |" << std::endl;
	ofs_running << " |                                            "
		"                        |" << std::endl;
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    assert(kpar>=1);
    assert(my_pool>=0);
    assert(my_rank>=0);
    assert(nbands>0);
    assert(nspin>0);
    assert(npol==1 || npol==2);

    const int nrxx = pw_wfc->nrxx;
    assert(nrxx>=0);

    for (int is = 0; is < nspin; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(chg.rho[is], nrxx);
    }

    const int ng_npol = pw_wfc->npwk_max * npol;

    ModuleBase::ComplexMatrix wfc_tmp(nbands, ng_npol);
    std::vector<std::complex<double>> rho_tmp(nrxx);

    // read occupation numbers
    ModuleBase::matrix wg_tmp(nkstot, nbands);
    if (my_rank == 0)
    {
        std::string filename = readin_dir + "eig.txt";
        std::ifstream ifs(filename);

		if(!ifs)
		{
            std::stringstream sss;
            sss << "Cannot find file " << filename;
			ModuleBase::WARNING_QUIT("ModuleIO::read_wf2rho_pw", sss.str());
		}
        else
        {
            ofs_running << " Find file containing weights of wave function: " << filename << std::endl;
        }

		std::string useless;
		getline(ifs, useless);
		getline(ifs, useless);
		for (int ik_tot = 0; ik_tot < nkstot; ++ik_tot)
		{
			ifs >> useless;
			getline(ifs, useless);
			for (int ib = 0; ib < nbands; ++ib)
			{
				ifs >> useless >> useless >> wg_tmp(ik_tot, ib);
			}
		}
	}

#ifdef __MPI
    MPI_Bcast(wg_tmp.c, nkstot * nbands, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    for (int ik = 0; ik < pw_wfc->nks; ++ik)
    {
        int is = 0;
        if (nspin == 2)
        {
            is = isk[ik];
        }
        const int ikstot = ik2iktot[ik];

        // mohan add 2025-05-17
        // .dat file
        const int out_type = 2;
        const bool out_app_flag = false;
        const bool gamma_only = false;
        const int istep = -1;

        std::string fn = filename_output(readin_dir,"wf","pw",ik,ik2iktot,nspin,nkstot,
                out_type,out_app_flag,gamma_only,istep);

        ofs_running << " Reading wave function from file: " << fn << std::endl;

		ModuleIO::read_wfc_pw(fn, pw_wfc, 
				rank_in_pool, nproc_in_pool, nbands, npol,
				ik, ikstot, nkstot, wfc_tmp);

        if (nspin == 4)
        {
            std::vector<std::complex<double>> rho_tmp2(nrxx);
            for (int ib = 0; ib < nbands; ++ib)
            {
                const std::complex<double>* wfc_ib = wfc_tmp.c + ib * ng_npol;
                const std::complex<double>* wfc_ib2 = wfc_tmp.c + ib * ng_npol + ng_npol / 2;
                pw_wfc->recip2real(wfc_ib, rho_tmp.data(), ik);
                pw_wfc->recip2real(wfc_ib2, rho_tmp2.data(), ik);
                const double w1 = wg_tmp(ikstot, ib) / pw_wfc->omega;

                if (w1 != 0.0)
                {
                    base_device::DEVICE_CPU* ctx = nullptr;
                    elecstate::elecstate_pw_op<double, base_device::DEVICE_CPU>()(ctx,
                                                                                  PARAM.globalv.domag,
                                                                                  PARAM.globalv.domag_z,
                                                                                  nrxx,
                                                                                  w1,
                                                                                  chg.rho,
                                                                                  rho_tmp.data(),
                                                                                  rho_tmp2.data());
                }
            }
        }
        else
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                const std::complex<double>* wfc_ib = wfc_tmp.c + ib * ng_npol;
                pw_wfc->recip2real(wfc_ib, rho_tmp.data(), ik);

                const double w1 = wg_tmp(ikstot, ib) / pw_wfc->omega;

                if (w1 != 0.0)
                {
					base_device::DEVICE_CPU* ctx = nullptr;
					elecstate::elecstate_pw_op<double, base_device::DEVICE_CPU>()(ctx, is, nrxx, 
							w1, chg.rho, rho_tmp.data());
                }
            }
        }
    }

#ifdef __MPI
    chg.init_chgmpi();
    for (int is = 0; is < nspin; ++is)
    {
        chg.reduce_diff_pools(chg.rho[is]);
    }
#endif

    // Since rho is calculated by psi^2, it is not symmetric. We need to rearrange it. 
    Symmetry_rho srho;
    for (int is = 0; is < nspin; is++)
    {
        srho.begin(is, chg, chg.rhopw, symm);
    }

    ModuleBase::timer::tick("ModuleIO", "read_wf2rho_pw");
}
