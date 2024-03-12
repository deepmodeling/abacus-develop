#include "write_wfc_pw.h"

#ifdef __MPI
#include "mpi.h"
#endif

#include "binstream.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include "module_base/tool_title.h"
#include "module_io/input.h"

void ModuleIO::write_wfc_pw(const std::string& fn,
                            const psi::Psi<std::complex<double>>& psi,
                            const K_Vectors& kv,
                            const ModulePW::PW_Basis_K* wfcpw)
{
    ModuleBase::TITLE("ModuleIO","write_wfc_pw");
    const int nkstot = kv.nkstot;
    const int nks = kv.nks;

    std::string * wfilename;
    wfilename = new std::string[nkstot];
    for(int ik = 0; ik < nkstot ; ++ik)
    { 
        std::stringstream wfss;
        if(INPUT.out_wfc_pw==1)
            wfss<<fn<<ik+1<<".txt";
        else if(INPUT.out_wfc_pw==2)
        {
            wfss<<fn<<ik+1<<".dat";
        }
        wfilename[ik]=wfss.str();
        if ( GlobalV::MY_RANK == 0 )
        {
            if(INPUT.out_wfc_pw==1)
            {
                std::ofstream ofs(wfss.str().c_str()); //clear all wavefunc files.
                ofs.close();
            }
            else if(INPUT.out_wfc_pw==2)
            {
                Binstream wfs(wfss.str(),"w");
                wfs.close();
            }
        }
    }
    //if(GlobalV::MY_RANK!=0) std::cout.clear();
    //std::cout<<"Hello"<<std::endl;
    //if(GlobalV::MY_RANK!=0) std::cout.setstate(ios::failbit);
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);

	
    // out put the wave functions in plane wave basis.
	for(int ip=0; ip<GlobalV::KPAR; ip++)
	{
        if( GlobalV::MY_POOL == ip )
		{
#endif
			for(int ik=0; ik<psi.get_nk(); ik++)
			{
                psi.fix_k(ik);
                int ikngtot=0; //ikngtot: the total number of plane waves of ikpoint
                int ikstot=0;//ikstot : the index within all k-points
#ifdef __MPI
                MPI_Allreduce(&kv.ngk[ik], &ikngtot, 1, MPI_INT, MPI_SUM, POOL_WORLD);

                // ikstot=GlobalC::Pkpoints.startk_pool[ip]+ik;
                // In the future, Pkpoints should be moved into Klist
                // To avoid GlobalC, we use getik_global instead
                ikstot = kv.getik_global(ik);
#else
                ikngtot = kv.ngk[ik];
                ikstot=ik;
#endif
#ifdef __MPI
				for( int id=0; id<GlobalV::NPROC_IN_POOL; id++)
				{
					MPI_Barrier(POOL_WORLD);
                    if (GlobalV::RANK_IN_POOL == id)
					{
#else
                    int id=0;               
#endif
                    if(INPUT.out_wfc_pw==1)
                    {
                        std::ofstream ofs2( wfilename[ikstot].c_str(),std::ios::app);
                        if(id==0)
                        {
                            ofs2<<std::setprecision(6);
                            ofs2<<std::setw(10)<<"Kpoint"<<std::setw(10)<<"nKpoint"<<std::setw(10)<<"kv.x"<<std::setw(10)
                            <<"kv.y"<<std::setw(10)<<"kv.z"<<std::setw(10)<<"weight"<<std::setw(10)
                            <<"ngtot"<<std::setw(10)<<"nband"<<std::setw(10)<<"ecut"<<std::setw(10)<<"lat0"<<std::setw(10)<<"2pi/lat0"<<std::endl;
                            ofs2 << std::setw(10) << ikstot + 1 << std::setw(10) << nkstot << std::setw(10)
                                 << kv.kvec_c[ik].x << std::setw(10) << kv.kvec_c[ik].y << std::setw(10)
                                 << kv.kvec_c[ik].z << std::setw(10) << kv.wk[ik] << std::setw(10) << ikngtot
                                 << std::setw(10) << GlobalV::NBANDS << std::setw(10) << INPUT.ecutwfc << std::setw(10)
                                 << wfcpw->lat0 << std::setw(10) << wfcpw->tpiba << std::endl;
                            ofs2<<"\n<Reciprocal Lattice Vector>"<<std::endl;
                            ofs2<<std::setw(10)<<wfcpw->G.e11<<std::setw(10)<<wfcpw->G.e12<<std::setw(10)<<wfcpw->G.e13<<std::endl;
                            ofs2<<std::setw(10)<<wfcpw->G.e21<<std::setw(10)<<wfcpw->G.e22<<std::setw(10)<<wfcpw->G.e23<<std::endl;
                            ofs2<<std::setw(10)<<wfcpw->G.e31<<std::setw(10)<<wfcpw->G.e32<<std::setw(10)<<wfcpw->G.e33<<std::endl;
                            ofs2<<"<Reciprocal Lattice Vector>\n"<<std::endl;
                            ofs2<<"<G vectors>"<<std::endl;
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
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            ofs2<<"<G vectors>\n"<<std::endl;
                        }
						ofs2.close();
                    }
                    else if(INPUT.out_wfc_pw==2)
                    {
                        Binstream wfs2( wfilename[ikstot],"a");
                        if(id==0)
                        {
                            wfs2 << int(72) << ikstot + 1 << nkstot << kv.kvec_c[ik].x << kv.kvec_c[ik].y
                                 << kv.kvec_c[ik].z << kv.wk[ik] << ikngtot << GlobalV::NBANDS << INPUT.ecutwfc
                                 << wfcpw->lat0 << wfcpw->tpiba << 72; // 4 int + 7 double is 72B
                            wfs2<<72<<wfcpw->G.e11<<wfcpw->G.e12<<wfcpw->G.e13
                                    <<wfcpw->G.e21<<wfcpw->G.e22<<wfcpw->G.e23
                                    <<wfcpw->G.e31<<wfcpw->G.e32<<wfcpw->G.e33<<72; //9 double is 72B
                        }
                        if(id==0)
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
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            wfs2 << ikngtot * 4 * 3;
                        } 
                        wfs2.close();
                    } 
#ifdef __MPI
					}
				}// end id
#endif
                for(int ib=0; ib<GlobalV::NBANDS; ib++)
				{
#ifdef __MPI
					for( int id=0; id<GlobalV::NPROC_IN_POOL; id++)
					{
						MPI_Barrier(POOL_WORLD); //qianrui add
                        if (GlobalV::RANK_IN_POOL == id)
						{
#else
                    int id=0;
#endif
                        if(INPUT.out_wfc_pw==1)
                        {
                            std::ofstream ofs2( wfilename[ikstot].c_str(),std::ios::app);
                            if(id==0)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2 << std::scientific;
							
                            for (int ig=0; ig<psi.get_current_nbas(); ig++)
							{
								if (ig%4==0&&(ig!=0||id!=0)) ofs2 << "\n";
								ofs2 << std::setw(15) << psi(ib, ig).real()
									<< std::setw(15) << psi(ib, ig).imag();
							} // end ig
                            if(id==GlobalV::NPROC_IN_POOL-1)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2.close();
                        }
                        else if(INPUT.out_wfc_pw==2)
                        {
                            Binstream wfs2(wfilename[ikstot],"a");
                            if(id==0) wfs2<<ikngtot*16;
                            for (int ig=0; ig<psi.get_current_nbas(); ig++)
							{
								wfs2 << psi(ib, ig).real() << psi(ib, ig).imag();
							}
                            if(id==GlobalV::NPROC_IN_POOL-1) wfs2<<ikngtot*16;
                            wfs2.close();
                        }
#ifdef __MPI                      
						}
					}// end id
#endif
				}// end ib
			}// end ik
#ifdef __MPI
		}
	}// end ip
#endif


    delete [] wfilename;
	return;
}
