#include "wf_io.h"
#include "../src_pw/global.h"
//#include "../../src_develop/src_wannier/eximport.h"

//qianrui add 2020-10-15
/*inline int getink(const int &ik,const int &ipool,const int &nktot,const int &kpar)
{
    int nkp=nktot/kpar;
    int rem=nktot%kpar;
    if(ipool<rem)
    {
        return ipool*nkp+ipool+ik;
    }
    else
    {
        return ipool*nkp+rem+ik;       
    }
    
}*/

void WF_io::write_wfc(const std::string &fn, const ModuleBase::ComplexMatrix *psi)
{
    if (GlobalV::test_wf) ModuleBase::TITLE("WF_io","write_wfc");

    std::ofstream ofs( fn.c_str() );

    //eximport exi;
    //exi.out_input(ofs);
    //exi.out_kpoints(ofs);
    //exi.out_igk(ofs);
    //exi.out_planewave(ofs);

    ofs << "\n<WAVEFUNC>";
    ofs << "\n" << GlobalV::NBANDS << " Number of bands." << std::endl;
    ofs << std::setprecision(6);
    for (int i=0; i<GlobalV::NBANDS; i++)
    {
        for (int ik=0; ik<GlobalC::kv.nks; ik++)
        {
            ofs << "\n" << ik;
            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
            {
                if (ig%4==0) ofs << "\n";
                ofs << std::setw(15) << psi[ik](i, ig).real()
                << std::setw(15) << psi[ik](i, ig).imag();
            }
            ofs << "\n";
        }
    }
    ofs << "\n<WAVEFUNC>";

    ofs.close();
    return;
}

void WF_io::write_wfc2(const std::string &fn, const ModuleBase::ComplexMatrix *psi, const ModuleBase::Vector3<double> *gkk)
{
    if (GlobalV::test_wf) ModuleBase::TITLE("WF_io","write_wfc2"); 

    std::string * wfilename;
    wfilename=new std::string[GlobalC::kv.nkstot];
    for(int ik=0;ik<GlobalC::kv.nkstot;ik++)
    {
        std::stringstream wfss;
        if(GlobalC::wf.out_wfc_pw==1)
            wfss<<fn<<ik+1<<".txt";
        else if(GlobalC::wf.out_wfc_pw==2)
        {
            wfss<<fn<<ik+1<<".dat";
        }
        wfilename[ik]=wfss.str();
        if ( GlobalV::MY_RANK == 0 )
        {
            if(GlobalC::wf.out_wfc_pw==1)
            {
                std::ofstream ofs(wfss.str().c_str()); //clear all wavefunc files.
                ofs.close();
            }
            else if(GlobalC::wf.out_wfc_pw==2)
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
			for(int ik=0; ik<GlobalC::kv.nks; ik++)
			{
                int ikngtot=0; //ikngtot: the total number of plane waves of ikpoint
                int ikstot=0;//ikstot : the index within all k-points
#ifdef __MPI
                MPI_Allreduce(&GlobalC::kv.ngk[ik],&ikngtot,1,MPI_INT,MPI_SUM,POOL_WORLD);
                ikstot=GlobalC::Pkpoints.startk_pool[ip]+ik;
#else
                ikngtot=GlobalC::kv.ngk[ik];
                ikstot=ik;
#endif
#ifdef __MPI
                //int ikstot=getink(ik,ip,GlobalC::kv.nkstot,GlobalV::KPAR);
				for( int id=0; id<GlobalV::NPROC_IN_POOL; id++)
				{
					MPI_Barrier(POOL_WORLD);
                    if (GlobalV::RANK_IN_POOL == id)
					{
#else
                    int id=0;               
#endif
                    if(GlobalC::wf.out_wfc_pw==1)
                    {
                        std::ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                        if(id==0)
                        {
                            ofs2<<std::setprecision(6);
                            ofs2<<std::setw(10)<<"Kpoint"<<std::setw(10)<<"nKpoint"<<std::setw(10)<<"GlobalC::kv.x"<<std::setw(10)
                            <<"GlobalC::kv.y"<<std::setw(10)<<"GlobalC::kv.z"<<std::setw(10)<<"weight"<<std::setw(10)
                            <<"ngtot"<<std::setw(10)<<"nband"<<std::setw(10)<<"ecut"<<std::setw(10)<<"lat0"<<std::setw(10)<<"2pi/lat0"<<std::endl;
                            ofs2<<std::setw(10)<<ikstot+1<<std::setw(10)<<GlobalC::kv.nkstot<<std::setw(10)<<GlobalC::kv.kvec_c[ik].x<<std::setw(10)
                            <<GlobalC::kv.kvec_c[ik].y<<std::setw(10)<<GlobalC::kv.kvec_c[ik].z<<std::setw(10)<<GlobalC::kv.wk[ik]<<std::setw(10)
                            <<ikngtot<<std::setw(10)<<GlobalV::NBANDS<<std::setw(10)<<GlobalC::pw.ecutwfc<<std::setw(10)<<GlobalC::ucell.lat0<<std::setw(10)<<GlobalC::ucell.tpiba<<std::endl;
                            ofs2<<"\n<Reciprocal Lattice Vector>"<<std::endl;
                            ofs2<<std::setw(10)<<GlobalC::ucell.G.e11<<std::setw(10)<<GlobalC::ucell.G.e12<<std::setw(10)<<GlobalC::ucell.G.e13<<std::endl;
                            ofs2<<std::setw(10)<<GlobalC::ucell.G.e21<<std::setw(10)<<GlobalC::ucell.G.e22<<std::setw(10)<<GlobalC::ucell.G.e23<<std::endl;
                            ofs2<<std::setw(10)<<GlobalC::ucell.G.e31<<std::setw(10)<<GlobalC::ucell.G.e32<<std::setw(10)<<GlobalC::ucell.G.e33<<std::endl;
                            ofs2<<"<Reciprocal Lattice Vector>\n"<<std::endl;
                            ofs2<<"<G vectors>"<<std::endl;
                        }
                        for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
					    {
                            ofs2<<std::setw(10)<<gkk[ig].x<<std::setw(10)<<gkk[ig].y<<std::setw(10)<<gkk[ig].z<<std::endl;
						}
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            ofs2<<"<G vectors>\n"<<std::endl;
                        }
						ofs2.close();
                    }
                    else if(GlobalC::wf.out_wfc_pw==2)
                    {
                        Rwstream wfs2( wfilename[ikstot],"a");
                        if(id==0)
                        {
                            wfs2<<int(72)<<ikstot+1<<GlobalC::kv.nkstot<<GlobalC::kv.kvec_c[ik].x
                                <<GlobalC::kv.kvec_c[ik].y<<GlobalC::kv.kvec_c[ik].z<<GlobalC::kv.wk[ik]
                                <<ikngtot<<GlobalV::NBANDS<<GlobalC::pw.ecutwfc<<GlobalC::ucell.lat0<<GlobalC::ucell.tpiba<<72; //4 int + 7 double is 72B
                            wfs2<<72<<GlobalC::ucell.G.e11<<GlobalC::ucell.G.e12<<GlobalC::ucell.G.e13
                                    <<GlobalC::ucell.G.e21<<GlobalC::ucell.G.e22<<GlobalC::ucell.G.e23
                                    <<GlobalC::ucell.G.e31<<GlobalC::ucell.G.e32<<GlobalC::ucell.G.e33<<72; //9 double is 72B
                        }
                        if(id==0)
                        {
                            wfs2<<ikngtot*8*3;
                        }
                         
                        for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
					    {
                            wfs2<<gkk[ig].x<<gkk[ig].y<<gkk[ig].z;
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
                        if(GlobalC::wf.out_wfc_pw==1)
                        {
                            std::ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                            if(id==0)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2 << scientific;
							
                            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
							{
								if (ig%4==0&&(ig!=0||id!=0)) ofs2 << "\n";
								ofs2 << std::setw(15) << psi[ik](ib, ig).real()
									<< std::setw(15) << psi[ik](ib, ig).imag();
							} // end ig
                            if(id==GlobalV::NPROC_IN_POOL-1)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2.close();
                        }
                        else if(GlobalC::wf.out_wfc_pw==2)
                        {
                            Rwstream wfs2(wfilename[ikstot],"a");
                            if(id==0) wfs2<<ikngtot*16;
                            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
							{
								wfs2 << psi[ik](ib, ig).real() << psi[ik](ib, ig).imag();
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
