#include "write_wfc_pw.h"

#ifdef __MPI
#include "mpi.h"
#endif

#include "binstream.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_global.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/filename.h"

void ModuleIO::write_wfc_pw(
        const int istep,
        const int iter,
        const int kpar,
        const int my_pool,
        const int my_rank,
        const int nbands,
        const int nspin,
        const int npol,
        const int rank_in_pool,
        const int nproc_in_pool,
        const int out_wfc_pw,
        const double& ecutwfc,
        const std::string& global_out_dir,
        const psi::Psi<std::complex<double>>& psi,
        const K_Vectors& kv,
        const ModulePW::PW_Basis_K* wfcpw,
        std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "write_wfc_pw");

	if(out_wfc_pw!=1 && out_wfc_pw!=2)
	{
		return;
	}

    const int nkstot = kv.get_nkstot();
    const int nks = kv.get_nks();

    assert(nkstot>0);
    assert(nks>0);

    bool out_app_flag = false; // need to modify later, mohan 2025-05-17
    bool gamma_only = false; // need to modify later, mohan 2025-05-17

    std::string* wfilename = new std::string[nks];

    for (int ip = 0; ip < kpar; ip++)
    {
        if (my_pool != ip) continue;

        for (int ik_local = 0; ik_local < kv.get_nks(); ik_local++)
        {
            std::string fn = filename_output(global_out_dir,"wf","pw",
                    ik_local,kv.ik2iktot,nspin,nkstot,
                    out_wfc_pw,out_app_flag,gamma_only,istep,iter);

            ofs_running << " Write G-space wave functions to file: "
                << fn << std::endl;

			wfilename[ik_local] = fn;

			if (rank_in_pool == 0)
			{
				if (out_wfc_pw == 1)
				{
					std::ofstream ofs(fn.c_str()); // clear all wavefunc files.
					ofs.close();
				}
				else if (out_wfc_pw == 2)
				{
					Binstream wfs(fn, "w");
					wfs.close();
				}
			}
		}
	}


#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif



#ifdef __MPI
    // out put the wave functions in plane wave basis.
    for (int ip = 0; ip < kpar; ip++)
    {
        if (my_pool == ip)
        {
#endif
            for (int ik = 0; ik < psi.get_nk(); ik++)
            {
                psi.fix_k(ik);
                int ikngtot = 0; // ikngtot: the total number of plane waves of ikpoint
                const int ng = kv.ngk[ik];
                const int ng_max = wfcpw->npwk_max;
                const int ikstot = kv.ik2iktot[ik];
#ifdef __MPI
                MPI_Allreduce(&kv.ngk[ik], &ikngtot, 1, MPI_INT, MPI_SUM, POOL_WORLD);
#else
                ikngtot = kv.ngk[ik];
#endif
                const int ikngtot_npol = ikngtot * npol;
#ifdef __MPI
                for (int id = 0; id < nproc_in_pool; id++)
                {
                    MPI_Barrier(POOL_WORLD);
                    if (rank_in_pool == id)
                    {
#else
        int id = 0;
#endif
                        if (out_wfc_pw == 1)
                        {
                            std::ofstream ofs2(wfilename[ik].c_str(), std::ios::app);
                            if (id == 0)
                            {
                                ofs2 << std::setprecision(6);
                                ofs2 << std::setw(10) << "Kpoint" << std::setw(10) << "nKpoint" << std::setw(10)
                                     << "kv.x" << std::setw(10) << "kv.y" << std::setw(10) << "kv.z" << std::setw(10)
                                     << "weight" << std::setw(10) << "ngtot" << std::setw(10) << "nband"
                                     << std::setw(10) << "ecut" << std::setw(10) << "lat0" << std::setw(10)
                                     << "2pi/lat0" << std::endl;
                                ofs2 << std::setw(10) << ikstot + 1 << std::setw(10) << nkstot << std::setw(10)
                                     << kv.kvec_c[ik].x << std::setw(10) << kv.kvec_c[ik].y << std::setw(10)
                                     << kv.kvec_c[ik].z << std::setw(10) << kv.wk[ik] << std::setw(10) << ikngtot
                                     << std::setw(10) << nbands << std::setw(10) << ecutwfc
                                     << std::setw(10) << wfcpw->lat0 << std::setw(10) << wfcpw->tpiba << std::endl;
                                ofs2 << "\n<Reciprocal Lattice Vector>" << std::endl;
                                ofs2 << std::setw(10) << wfcpw->G.e11 << std::setw(10) << wfcpw->G.e12 << std::setw(10)
                                     << wfcpw->G.e13 << std::endl;
                                ofs2 << std::setw(10) << wfcpw->G.e21 << std::setw(10) << wfcpw->G.e22 << std::setw(10)
                                     << wfcpw->G.e23 << std::endl;
                                ofs2 << std::setw(10) << wfcpw->G.e31 << std::setw(10) << wfcpw->G.e32 << std::setw(10)
                                     << wfcpw->G.e33 << std::endl;
                                ofs2 << "<Reciprocal Lattice Vector>\n" << std::endl;
                                ofs2 << "<G vectors>" << std::endl;
                            }
                            for (int igl = 0; igl < wfcpw->npwk[ik]; ++igl)
                            {
                                int isz = wfcpw->igl2isz_k[ik * wfcpw->npwk_max + igl];
                                int iz = isz % wfcpw->nz;
                                int is = isz / wfcpw->nz;
                                int ixy = wfcpw->is2fftixy[is];
                                int ix = ixy / wfcpw->fftny;
                                int iy = ixy % wfcpw->fftny;

                                ofs2 << std::setw(10) << ix << std::setw(10) << iy << std::setw(10) << iz << std::endl;
                            }
                            if (id == nproc_in_pool - 1)
                            {
                                ofs2 << "<G vectors>\n" << std::endl;
                            }
                            ofs2.close();
                        }
                        else if (out_wfc_pw == 2)
                        {
                            Binstream wfs2(wfilename[ik], "a");
                            if (id == 0)
                            {
                                wfs2 << int(72) << ikstot + 1 << nkstot << kv.kvec_c[ik].x << kv.kvec_c[ik].y
                                     << kv.kvec_c[ik].z << kv.wk[ik] << ikngtot << nbands << ecutwfc
                                     << wfcpw->lat0 << wfcpw->tpiba << 72; // 4 int + 7 double is 72B
                                wfs2 << 72 << wfcpw->G.e11 << wfcpw->G.e12 << wfcpw->G.e13 << wfcpw->G.e21
                                     << wfcpw->G.e22 << wfcpw->G.e23 << wfcpw->G.e31 << wfcpw->G.e32 << wfcpw->G.e33
                                     << 72; // 9 double is 72B
                            }
                            if (id == 0)
                            {
                                wfs2 << ikngtot * 4 * 3;
                            }
                            for (int igl = 0; igl < wfcpw->npwk[ik]; ++igl)
                            {
                                int isz = wfcpw->igl2isz_k[ik * wfcpw->npwk_max + igl];
                                int iz = isz % wfcpw->nz;
                                int is = isz / wfcpw->nz;
                                int ixy = wfcpw->is2fftixy[is];
                                int ix = ixy / wfcpw->fftny;
                                int iy = ixy % wfcpw->fftny;

                                wfs2 << ix << iy << iz;
                            }
                            if (id == nproc_in_pool - 1)
                            {
                                wfs2 << ikngtot * 4 * 3;
                            }
                            wfs2.close();
                        }
#ifdef __MPI
                    }
                } // end id
#endif
                for (int ib = 0; ib < nbands; ib++)
                {
#ifdef __MPI
                    for (int id = 0; id < nproc_in_pool; id++)
                    {
                        MPI_Barrier(POOL_WORLD); // qianrui add
                        if (rank_in_pool == id)
                        {
#else
                        int id = 0;
#endif
                            if (out_wfc_pw == 1)
                            {
                                std::ofstream ofs2(wfilename[ik].c_str(), std::ios::app);
                                if (id == 0)
                                {
                                    ofs2 << "\n< Band " << ib + 1 << " >" << std::endl;
                                }
                                ofs2 << std::scientific;
                                for (int ig = 0; ig < ng; ig++)
                                {
                                    if (ig % 4 == 0 && (ig != 0 || id != 0))
                                    {
                                        ofs2 << "\n";
                                    }
                                    ofs2 << std::setw(15) << psi(ib, ig).real() << std::setw(15) << psi(ib, ig).imag();
                                } // end ig
                                if (id == nproc_in_pool - 1 && npol == 1)
                                {
                                    ofs2 << "\n< Band " << ib + 1 << " >" << std::endl;
                                }
                                ofs2.close();
                            }
                            else if (out_wfc_pw == 2)
                            {
                                Binstream wfs2(wfilename[ik], "a");
                                if (id == 0)
                                {
                                    wfs2 << ikngtot_npol * 16;
                                }
                                for (int ig = 0; ig < ng; ig++)
                                {
                                    wfs2 << psi(ib, ig).real() << psi(ib, ig).imag();
                                }
                                if (id == nproc_in_pool - 1 && npol == 1)
                                {
                                    wfs2 << ikngtot_npol * 16;
                                }
                                wfs2.close();
                            }
#ifdef __MPI
                        }
                    } // end id
#endif
                    if (npol > 1)
                    {
#ifdef __MPI
                        for (int id = 0; id < nproc_in_pool; id++)
                        {
                            MPI_Barrier(POOL_WORLD); // qianrui add
                            if (rank_in_pool == id)
                            {
#else
                        int id = 0;
#endif
                                if (out_wfc_pw == 1)
                                {
                                    std::ofstream ofs2(wfilename[ik].c_str(), std::ios::app);

                                    ofs2 << std::scientific;
                                    for (int ig = 0; ig < ng; ig++)
                                    {
                                        if (ig % 4 == 0 && (ig != 0 || id != 0))
                                        {
                                            ofs2 << "\n";
                                        }
                                        ofs2 << std::setw(15) << psi(ib, ig + ng_max).real() << std::setw(15)
                                             << psi(ib, ig + ng_max).imag();
                                    } // end ig
                                    if (id == nproc_in_pool - 1)
                                    {
                                        ofs2 << "\n< Band " << ib + 1 << " >" << std::endl;
                                    }
                                    ofs2.close();
                                }
                                else if (out_wfc_pw == 2)
                                {
                                    Binstream wfs2(wfilename[ik], "a");
                                    for (int ig = 0; ig < ng; ig++)
                                    {
                                        wfs2 << psi(ib, ig + ng_max).real() << psi(ib, ig + ng_max).imag();
                                    }
                                    if (id == nproc_in_pool - 1)
                                    {
                                        wfs2 << ikngtot_npol * 16;
                                    }
                                    wfs2.close();
                                }
#ifdef __MPI
                            } // end if rank_in_pool
                        } // end id
#endif
                    } // end if npol>1

                } // end ib
            }     // end ik
#ifdef __MPI
        } // end if my_pool
    }     // end ipool
#endif

    delete[] wfilename;
    return;
}
