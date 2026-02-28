#include "write_bands.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_base/formatter.h"
#include "source_io/module_parameter/parameter.h" // mohan add 20250911

#ifdef __MPI
#include <mpi.h>
#endif


void ModuleIO::write_bands(const Input_para& inp,
    const ModuleBase::matrix& ekb,
    const K_Vectors& kv)
{
    // write band information to band.txt
    if (inp.out_band[0])
    {
        const int nspin0 = (inp.nspin == 2) ? 2 : 1;
        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "band";

            if(nspin0==1)
            {
                // do nothing
            }
            else if(nspin0==2)
            {
                ss << "s" << is + 1;
            }

            ss << ".txt";

            const double eshift = 0.0;
            nscf_bands(is, ss.str(), inp.nbands, eshift, 
                       inp.out_band[1], // precision
                       ekb, kv);
        }
    }
}


void ModuleIO::nscf_bands(
    const int &is,
    const std::string &eig_file, 
    const int &nband,
    const double &fermie,
    const int &precision,
    const ModuleBase::matrix& ekb,
    const K_Vectors& kv)
{
    ModuleBase::TITLE("ModuleIO","nscf_bands");
    ModuleBase::timer::tick("ModuleIO", "nscf_bands");

    assert(precision>0);

/*
	GlobalV::ofs_running << "\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " |             #Print out the eigenvalues for each spin#              |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n";
*/

    GlobalV::ofs_running << " Write eigenvalues to file: " << eig_file << std::endl;

    // number of k points without spin; 
    // nspin = 1,2, nkstot = nkstot_np * nspin; 
    // nspin = 4, nkstot = nkstot_np
    const int nkstot_np = kv.para_k.nkstot_np;
    const int nks_np = kv.para_k.nks_np;


#ifdef __MPI
    if(GlobalV::MY_RANK==0)
	{
		if(PARAM.inp.out_app_flag==true)
		{
            // do nothing
		}
		else
		{
			std::ofstream ofs(eig_file.c_str());
			ofs.close();
		}
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<double> klength;
    klength.resize(nkstot_np);
    klength[0] = 0.0;
    std::vector<ModuleBase::Vector3<double>> kvec_c_global;
    kv.para_k.gatherkvec(kv.kvec_c, kvec_c_global);

    for(int ik=0; ik<nkstot_np; ik++)
    {
        if (ik>0)
        {
            auto delta=kvec_c_global[ik]-kvec_c_global[ik-1];
            klength[ik] = klength[ik-1];
            klength[ik] += (kv.kl_segids[ik] == kv.kl_segids[ik-1]) ? delta.norm() : 0.0;
        }
        //! first find if present kpoint in present pool
        if ( GlobalV::MY_POOL == kv.para_k.whichpool[ik] )
        {
            //! then get the local kpoint index, which starts definitly from 0
            const int ik_now = ik - kv.para_k.startk_pool[GlobalV::MY_POOL];
            //! if present kpoint corresponds the spin of the present one
            assert( kv.isk[ik_now+is*nks_np] == is );
            if ( GlobalV::RANK_IN_POOL == 0)
            {
                std::ofstream ofs(eig_file.c_str(), std::ios::app);
                ofs << FmtCore::format("%4d", ik+1);
                int width = precision + 4;
                std::string fmtstr = " %." + std::to_string(precision) + "f";
                ofs << FmtCore::format(fmtstr.c_str(), klength[ik]);
                for(int ib = 0; ib < nband; ib++)
                {
                    ofs << FmtCore::format(fmtstr.c_str(), (ekb(ik_now+is*nks_np, ib)-fermie) * ModuleBase::Ry_to_eV);
                }
                ofs << std::endl;
                ofs.close();    
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
#else
    std::vector<double> klength;
    klength.resize(nkstot_np);
    klength[0] = 0.0;

    std::ofstream ofs;

    if(GlobalV::MY_RANK==0)
	{
		if(PARAM.inp.out_app_flag==true)
		{
            ofs.open(eig_file.c_str(), std::ios::app);
		}
		else
		{
			ofs.open(eig_file.c_str());
		}
    }

    for(int ik=0;ik<nkstot_np;ik++)
    {
        if (ik>0)
        {
            auto delta=kv.kvec_c[ik]-kv.kvec_c[ik-1];
            klength[ik] = klength[ik-1];
            klength[ik] += (kv.kl_segids[ik] == kv.kl_segids[ik-1]) ? delta.norm() : 0.0;
        }
        if( kv.isk[ik] == is)
        {
            ofs << FmtCore::format("%4d", ik+1);
            int width = precision + 4;
            std::string fmtstr = " %." + std::to_string(precision) + "f";
            ofs << FmtCore::format(fmtstr.c_str(), klength[ik]);
            for(int ibnd = 0; ibnd < nband; ibnd++)
            {
                ofs << FmtCore::format(fmtstr.c_str(), (ekb(ik, ibnd)-fermie) * ModuleBase::Ry_to_eV);
            }
            ofs << std::endl;
        }
    }

    ofs.close();
#endif

    ModuleBase::timer::tick("ModuleIO", "nscf_band");
    return;
}
