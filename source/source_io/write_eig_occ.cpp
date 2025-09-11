#include "write_eig_occ.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/parallel_comm.h" // use POOL_WORLD

#ifdef __MPI
#include <mpi.h> // use MPI_Barrier
#endif

void ModuleIO::write_eig_iter(const ModuleBase::matrix &ekb,const ModuleBase::matrix &wg, const K_Vectors& kv)
{
    ModuleBase::TITLE("ModuleIO","write_eig_iter");
	ModuleBase::timer::tick("ModuleIO", "write_eig_iter");

	GlobalV::ofs_running << "\n PRINT #EIGENVALUES# AND #OCCUPATIONS#" << std::endl;

    const int nspin = PARAM.inp.nspin;
    const int nks = kv.get_nks();
	const int nkstot = kv.get_nkstot();

    std::vector<int> ngk_tot = kv.ngk;

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ngk_tot.data(), nks, MPI_INT, MPI_SUM, POOL_WORLD);
#endif

    const int nk_fac = nspin == 2 ? 2 : 1;
    const int nks_np = nks / nk_fac;
    const int nkstot_np = nkstot / nk_fac;
    const int kpar = GlobalV::KPAR;

    for (int is = 0; is < nk_fac; ++is)
    {
        for (int ip = 0; ip < kpar; ++ip)
        {
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            bool ip_flag = PARAM.inp.out_alllog || (GlobalV::RANK_IN_POOL == 0 && GlobalV::MY_BNDGROUP == 0);

            if (GlobalV::MY_POOL == ip && ip_flag)
            {
                GlobalV::ofs_running << std::setprecision(8);
//              ofs_eig << std::setiosflags(std::ios::showpoint);

                const int start_ik = nks_np * is;
                const int end_ik = nks_np * (is + 1);
                for (int ik = start_ik; ik < end_ik; ++ik)
                {
                    GlobalV::ofs_running << " spin=" << is+1 << " k-point="
                            << kv.ik2iktot[ik] + 1 - is * nkstot_np << std::endl;

                    GlobalV::ofs_running << std::setw(8) << "Index"
                    << std::setw(18) << "Eigenvalues(eV)"
                    << std::setw(18) << "Occupations" << std::endl;

                    for (int ib = 0; ib < ekb.nc; ib++)
                    {
                        GlobalV::ofs_running << std::setw(8) << ib + 1 
                                << std::setw(18) << ekb(ik, ib) * ModuleBase::Ry_to_eV
                                << std::setw(18) << wg(ik, ib) << std::endl;
                    }
                    GlobalV::ofs_running << std::endl;
                }
            }

	    
// =============================================================
// MPI communication: 
// RANK 0 collect and print out #EIGENVALUES# AND #OCCUPATIONS# for all k-points
// =============================================================
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);

            if (kpar > 1 && ip > 0 && PARAM.inp.out_alllog == 0)
            {
                        
                if (GlobalV::MY_RANK == 0 )
                {
                    
                    // for the current spin channel [is] and pool [ip] 
                    // MPI_Recv the size of matrix, ik2iktot, ekb and wg from RANK_IN_POOL=0
                    MPI_Status status;
                    int recv_nks_np, recv_nbands;
                    MPI_Recv(&recv_nks_np, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    int source_rank = status.MPI_SOURCE;
                    MPI_Recv(&recv_nbands, 1, MPI_INT, source_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int* recv_ik2iktot = new int[recv_nks_np]; 
                    MPI_Recv(recv_ik2iktot, recv_nks_np, MPI_INT, source_rank, 
                             2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    ModuleBase::matrix recv_ekb(recv_nks_np, recv_nbands);
                    MPI_Recv(recv_ekb.c, recv_nks_np * recv_nbands, MPI_DOUBLE, 
                             source_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    ModuleBase::matrix recv_wg(recv_nks_np, recv_nbands);
                    MPI_Recv(recv_wg.c, recv_nks_np * recv_nbands, MPI_DOUBLE, source_rank, 4, 
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                    
                    // print EIGENVALUES and OCCUPATIONS of received k-points
                    for (int ik = 0; ik < recv_nks_np; ++ik)
                    {
                        GlobalV::ofs_running << " spin=" << is+1 << " k-point="
                                << recv_ik2iktot[ik] + 1 - is * nkstot_np  << std::endl;

                        GlobalV::ofs_running << std::setw(8) << "Index"
                        << std::setw(18) << "Eigenvalues(eV)"
                        << std::setw(18) << "Occupations" << std::endl;

                        for (int ib = 0; ib < recv_nbands; ib++)
                        {
                            GlobalV::ofs_running << std::setw(8) << ib + 1 
                                    << std::setw(18) << recv_ekb(ik, ib) * ModuleBase::Ry_to_eV
                                    << std::setw(18) << recv_wg(ik, ib) << std::endl;
                        }
                        GlobalV::ofs_running << std::endl;
                    }
                            
                    delete[] recv_ik2iktot;
                }
                else if (GlobalV::MY_POOL == ip  && GlobalV::RANK_IN_POOL == 0 && GlobalV::MY_BNDGROUP == 0)
                {
                    // for the current spin channel [is] and pool [ip] 
                    // MPI_Send the size of matrix, ik2iktot, ekb and wg to RANK=0
                    const int send_nks_np = nks_np;          
                    const int send_nbands = ekb.nc;          
                    const int is_offset = is * nks_np;     
                    int* send_ik2iktot = new int[send_nks_np];
                    for (int ik = 0; ik < send_nks_np; ++ik)
                    {     
                        send_ik2iktot[ik] = kv.ik2iktot[is_offset + ik];    
                    }
                    MPI_Send(&send_nks_np, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);    
                    MPI_Send(&send_nbands, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    
                    MPI_Send(send_ik2iktot, send_nks_np, MPI_INT, 0, 2, MPI_COMM_WORLD); 
                    MPI_Send(ekb.c + is_offset * send_nbands,   
                            send_nks_np * send_nbands,           
                            MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);   
                    MPI_Send(wg.c + is_offset * send_nbands,    
                            send_nks_np * send_nbands,
                            MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);   

                    delete[] send_ik2iktot;
                }

            MPI_Barrier(MPI_COMM_WORLD);
            }
// =============================================================
// End of MPI Communication
// =============================================================
#endif
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    
	ModuleBase::timer::tick("ModuleIO", "write_eig_iter");
}

void ModuleIO::write_eig_file(const ModuleBase::matrix &ekb,
		const ModuleBase::matrix &wg, 
		const K_Vectors& kv,
		const int istep)
{
	ModuleBase::TITLE("ModuleIO","write_eig_file");
	ModuleBase::timer::tick("ModuleIO", "write_eig_file");

/*
	GlobalV::ofs_running << "\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " |            #Print out the eigenvalues and occupations#             |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n";
*/

    const int nspin = PARAM.inp.nspin;
    const int nks = kv.get_nks();
	const int nkstot = kv.get_nkstot();

    bool wrong = false;

	for (int ik = 0; ik < nks; ++ik)
	{
		for (int ib = 0; ib < ekb.nc; ++ib)
		{
			if (std::abs(ekb(ik, ib)) > 1.0e10)
			{
				GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1
					<< " " << ekb(ik, ib) << " Ry" << std::endl;
				wrong = true;
			}
		}
	}

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &wrong, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#endif
    if (wrong)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_eig_file", "Eigenvalues are too large!");
    }

    std::vector<int> ngk_tot = kv.ngk;

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ngk_tot.data(), nks, MPI_INT, MPI_SUM, POOL_WORLD);
#endif    

    // file name to store eigenvalues
    std::string filename = PARAM.globalv.global_out_dir + "eig_occ.txt";

    GlobalV::ofs_running << " Write eigenvalues and occupations to file: " << filename << std::endl;

    if (GlobalV::MY_RANK == 0)
    {
        std::ofstream ofs_eig0;

		if(PARAM.inp.out_app_flag==true)
		{
			ofs_eig0.open(filename.c_str(), std::ios::app);
		}
		else
		{
			ofs_eig0.open(filename.c_str());
		}
       
        ofs_eig0 << istep+1 << "     # ionic step" << std::endl;
        ofs_eig0 << " Electronic state energy (eV) and occupations" << std::endl;
        ofs_eig0 << " Spin number " << nspin << std::endl;
        ofs_eig0.close();
    }

    const int nk_fac = (nspin == 2) ? 2 : 1;
    const int nks_np = nks / nk_fac;
    const int nkstot_np = nkstot / nk_fac;
    const int kpar = GlobalV::KPAR;

    for (int is = 0; is < nk_fac; ++is)
    {
        for (int ip = 0; ip < kpar; ++ip)
        {
#ifdef __MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            bool ip_flag = PARAM.inp.out_alllog || (GlobalV::RANK_IN_POOL == 0 && GlobalV::MY_BNDGROUP == 0);

            if (GlobalV::MY_POOL == ip && ip_flag)
            {
                std::ofstream ofs_eig(filename.c_str(), std::ios::app);
                ofs_eig << std::setprecision(8);
                ofs_eig << std::setiosflags(std::ios::showpoint);

                const int start_ik = nks_np * is;
                const int end_ik = nks_np * (is + 1);
                for (int ik = start_ik; ik < end_ik; ++ik)
                {
                    ofs_eig << " spin=" << is+1 << " k-point="
                            << kv.ik2iktot[ik] + 1 - is * nkstot_np << "/" << nkstot_np
                            << " Cartesian=" << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y
                            << " " << kv.kvec_c[ik].z << " (" << ngk_tot[ik] << " plane wave)" << std::endl;

                    ofs_eig << std::setprecision(16);
                    ofs_eig << std::setiosflags(std::ios::showpoint);
                    for (int ib = 0; ib < ekb.nc; ib++)
                    {
                        ofs_eig << " " << ib + 1 << " " << ekb(ik, ib) * ModuleBase::Ry_to_eV
                                << " " << wg(ik, ib) << std::endl;
                    }
                    ofs_eig << std::endl;
                }

                ofs_eig.close();
            }
        }
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

	ModuleBase::timer::tick("ModuleIO", "write_eig_file");
	return;
}
