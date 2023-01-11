#include "wf_io.h"
#include "rwstream.h"
#include "module_base/tool_title.h"
#include "module_base/global_variable.h"
#include "input.h"
#include "../src_parallel/parallel_global.h"
#ifdef __MPI
#include "mpi.h"
#endif


inline int getink(const int &ik,const int &rankp,const int &nktot,const int &kpar)
{
    int nkp = nktot / kpar;
    int rem = nktot % kpar;
    if(rankp < rem)
    {
        return rankp*nkp + rankp + ik;
    }
    else
    {
        return rankp*nkp + rem + ik;       
    }
    
}

void WF_io::write_wfc(  const std::string &fn, const psi::Psi<std::complex<double>> &psi,
                        const K_Vectors* p_kv, const ModulePW::PW_Basis_K *wfc_basis)
{
    ModuleBase::TITLE("WF_io","write_wfc"); 
    const int nkstot = p_kv->nkstot;
    const int nks = p_kv->nks;

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
                Rwstream wfs(wfss.str(),"w");
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
                MPI_Allreduce(&p_kv->ngk[ik],&ikngtot,1,MPI_INT,MPI_SUM,POOL_WORLD);
                
                // ikstot=GlobalC::Pkpoints.startk_pool[ip]+ik;
                // In the future, Pkpoints should be moved into Klist
                // To avoid GlobalC, we use getink instead
                ikstot = getink(ik, ip, nkstot, GlobalV::KPAR);
#else
                ikngtot=p_kv->ngk[ik];
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
                        std::ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                        if(id==0)
                        {
                            ofs2<<std::setprecision(6);
                            ofs2<<std::setw(10)<<"Kpoint"<<std::setw(10)<<"nKpoint"<<std::setw(10)<<"kv.x"<<std::setw(10)
                            <<"kv.y"<<std::setw(10)<<"kv.z"<<std::setw(10)<<"weight"<<std::setw(10)
                            <<"ngtot"<<std::setw(10)<<"nband"<<std::setw(10)<<"ecut"<<std::setw(10)<<"lat0"<<std::setw(10)<<"2pi/lat0"<<std::endl;
                            ofs2<<std::setw(10)<<ikstot+1<<std::setw(10)<<nkstot<<std::setw(10)<<p_kv->kvec_c[ik].x<<std::setw(10)
                            <<p_kv->kvec_c[ik].y<<std::setw(10)<<p_kv->kvec_c[ik].z<<std::setw(10)<<p_kv->wk[ik]<<std::setw(10)
                            <<ikngtot<<std::setw(10)<<GlobalV::NBANDS<<std::setw(10)<<INPUT.ecutwfc<<std::setw(10)<<wfc_basis->lat0<<std::setw(10)<<wfc_basis->tpiba<<std::endl;
                            ofs2<<"\n<Reciprocal Lattice Vector>"<<std::endl;
                            ofs2<<std::setw(10)<<wfc_basis->G.e11<<std::setw(10)<<wfc_basis->G.e12<<std::setw(10)<<wfc_basis->G.e13<<std::endl;
                            ofs2<<std::setw(10)<<wfc_basis->G.e21<<std::setw(10)<<wfc_basis->G.e22<<std::setw(10)<<wfc_basis->G.e23<<std::endl;
                            ofs2<<std::setw(10)<<wfc_basis->G.e31<<std::setw(10)<<wfc_basis->G.e32<<std::setw(10)<<wfc_basis->G.e33<<std::endl;
                            ofs2<<"<Reciprocal Lattice Vector>\n"<<std::endl;
                            ofs2<<"<G vectors>"<<std::endl;
                        }
                        for (int ig=0; ig<p_kv->ngk[ik]; ig++)
					    {
                            ofs2<<std::setw(10)<<wfc_basis->getgcar(ik,ig).x<<std::setw(10)<<wfc_basis->getgcar(ik,ig).y<<std::setw(10)<<wfc_basis->getgcar(ik,ig).z<<std::endl;
						}
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            ofs2<<"<G vectors>\n"<<std::endl;
                        }
						ofs2.close();
                    }
                    else if(INPUT.out_wfc_pw==2)
                    {
                        Rwstream wfs2( wfilename[ikstot],"a");
                        if(id==0)
                        {
                            wfs2<<int(72)<<ikstot+1<<nkstot<<p_kv->kvec_c[ik].x
                                <<p_kv->kvec_c[ik].y<<p_kv->kvec_c[ik].z<<p_kv->wk[ik]
                                <<ikngtot<<GlobalV::NBANDS<<INPUT.ecutwfc<<wfc_basis->lat0<<wfc_basis->tpiba<<72; //4 int + 7 double is 72B
                            wfs2<<72<<wfc_basis->G.e11<<wfc_basis->G.e12<<wfc_basis->G.e13
                                    <<wfc_basis->G.e21<<wfc_basis->G.e22<<wfc_basis->G.e23
                                    <<wfc_basis->G.e31<<wfc_basis->G.e32<<wfc_basis->G.e33<<72; //9 double is 72B
                        }
                        if(id==0)
                        {
                            wfs2<<ikngtot*8*3;
                        }
                         
                        for (int ig=0; ig<p_kv->ngk[ik]; ig++)
					    {
                            wfs2<<wfc_basis->getgcar(ik,ig).x<<wfc_basis->getgcar(ik,ig).y<<wfc_basis->getgcar(ik,ig).z;
						}
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            wfs2<<ikngtot*8*3;
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
                            std::ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                            if(id==0)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2 << scientific;
							
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
                            Rwstream wfs2(wfilename[ikstot],"a");
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
