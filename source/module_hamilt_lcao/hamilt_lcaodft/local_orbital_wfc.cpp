#include "local_orbital_wfc.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_wfc_nao.h"
#include "module_io/read_wfc_nao.h"
#include "module_base/parallel_common.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

Local_Orbital_wfc::Local_Orbital_wfc()
{
	allocate_flag = false;
	wfck_flag = false;	
	complex_flag = false;
    nks = 0;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{

	// used for k-points.
	if(this->complex_flag)
	{
        delete[] this->wfc_k_grid2;
    }
    if(this->wfck_flag)
    {
		for(int i=0; i<nks; i++)
		{
			delete[] this->wfc_k_grid[i];
		}
		delete[] this->wfc_k_grid;
	}
}

void Local_Orbital_wfc::allocate_k(const int& lgd,
    psi::Psi<std::complex<double>>* psi,
    elecstate::ElecState* pelec,
    const int& nks,
    const int& nkstot,
    const std::vector<ModuleBase::Vector3<double>>& kvec_c)
{
    this->nks = nks;

    ModuleBase::TITLE("Local_Orbital_wfc", "allocate_k");
    if(GlobalV::NLOCAL < GlobalV::NBANDS)
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}

	// mohan add the flag 2011-03-02
	// allocate the first part (only once!).
	if(this->wfck_flag == false)
	{
		this->wfc_k_grid = new std::complex<double>**[nks];
		for(int ik=0; ik<nks; ik++)
		{
			this->wfc_k_grid[ik] = new std::complex<double>*[GlobalV::NBANDS];
		}
		this->wfck_flag = true;
	}
	
	if(this->complex_flag)
	{
		delete[] this->wfc_k_grid2;
		this->complex_flag = false;
	}
	// allocate the second part.
	//if(lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
	//if(lgd != 0 && this->complex_flag == false)
	if(lgd != 0)
	{
		//std::cout<<"lgd="<<lgd<<" ; GlobalV::NLOCAL="<<GlobalV::NLOCAL<<std::endl; //delete 2015-09-06, xiaohui
		const int page=GlobalV::NBANDS*lgd;
		this->wfc_k_grid2=new std::complex<double> [nks*page];
		ModuleBase::GlobalFunc::ZEROS(wfc_k_grid2, nks*page);
		for(int ik=0; ik<nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				this->wfc_k_grid[ik][ib] = &wfc_k_grid2[ik*page+ib*lgd];
				//std::cout<<"ik="<<ik<<" ib="<<ib<<std::endl<<"wfc_k_grid address: "<<wfc_k_grid[ik][ib]<<" wfc_k_grid2 address: "<<&wfc_k_grid2[ik*page+ib*lgd]<<std::endl;
			}
			//std::cout<<"set wfc_k_grid pointer success, ik: "<<ik<<std::endl;
			ModuleBase::Memory::record("LOWF::wfc_k_grid", sizeof(std::complex<double>) * GlobalV::NBANDS*GlobalV::NLOCAL);
			//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MemoryForWaveFunctions (MB)",mem);
			//std::cout<<"wfc_k_grid["<<ik<<"] use "<<mem<<" MB"<<std::endl;
			this->complex_flag = true;
		}
	}

    if (INPUT.init_wfc == "atomic")
    {
    }
    else if (INPUT.init_wfc == "file")
    {
        int error;
        std::cout << " Read in wave functions files: " << nkstot << std::endl;
        if(psi == nullptr)
        {
            ModuleBase::WARNING_QUIT("allocate_k","psi should be allocated first!");
        }
        else
        {
            psi->resize(nkstot, this->ParaV->ncol_bands, this->ParaV->nrow);
        }
		for(int ik=0; ik<nkstot; ++ik)
		{
            GlobalV::ofs_running << " Read in wave functions " << ik + 1 << std::endl;
            std::complex<double>** ctot;
            error = ModuleIO::read_wfc_nao_complex(ctot, ik, kvec_c[ik], this->ParaV, psi, pelec);
#ifdef __MPI
            Parallel_Common::bcast_int(error);
#endif
            GlobalV::ofs_running << " Error=" << error << std::endl;
            if(error==1)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc","Can't find the wave function file: LOWF.dat");
            }
            else if(error==2)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, band number doesn't match");
            }
            else if(error==3)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In wave function file, nlocal doesn't match");
            }
            else if(error==4)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc","In k-dependent wave function file, k point is not correct");
            }
        }
    }
    else
    {
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc","check the parameter: init_wfc");
	}

	return;
}
int Local_Orbital_wfc::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localindex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localindex%nblk;
    return gIndex;
}


int Local_Orbital_wfc::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalindex%(nblk*nprocs))/nblk);
    return int(globalindex/(nblk*nprocs))*nblk+globalindex%nblk;
}

#ifdef __MPI
void Local_Orbital_wfc::wfc_2d_to_grid(const int istep,
                                       int out_wfc_lcao,
                                       const double* wfc_2d,
                                       double** wfc_grid,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick(" Local_Orbital_wfc","wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int  myid;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info;
    
    //calculate maxnloc for bcasting 2d-wfc
    long maxnloc; // maximum number of elements in local matrix
    info=MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info=MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    double *work=new double[maxnloc]; // work/buffer matrix

    double** ctot;
    if (out_wfc_lcao && myid == 0)
    {
        ctot = new double* [GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new double[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
        }
        ModuleBase::Memory::record("LOWF::ctot", sizeof(double) * GlobalV::NBANDS * GlobalV::NLOCAL);
    }

    int naroc[2]; // maximum number of row or column
    for(int iprow=0; iprow<pv->dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            info=MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if(myid==src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work, inc);
                naroc[0]=pv->nrow;
                naroc[1]=pv->ncol_bands;
            }
            info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

            if (out_wfc_lcao)
                info = this->set_wfc_grid(naroc, pv->nb,
                    pv->dim0, pv->dim1, iprow, ipcol,
                    work, wfc_grid, myid, ctot);
            else
                info = this->set_wfc_grid(naroc, pv->nb,
                        pv->dim0, pv->dim1, iprow, ipcol,
                        work, wfc_grid);

        }//loop ipcol
    }//loop iprow
    if(out_wfc_lcao && myid == 0)
    {
        std::stringstream ss;
        if (GlobalV::out_app_flag)
        {
            ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".dat";
        }
        else
        {
            ss << GlobalV::global_out_dir << istep << "_"
               << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".dat";
        }
        ModuleIO::write_wfc_nao(ss.str(), ctot, ekb, wg);
        for (int i = 0; i < GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }
    delete[] work;
    ModuleBase::timer::tick(" Local_Orbital_wfc","wfc_2d_to_grid");
}

void Local_Orbital_wfc::wfc_2d_to_grid(const int istep,
                                       int out_wfc_lcao,
                                       const std::complex<double>* wfc_2d,
                                       std::complex<double>** wfc_grid,
                                       int ik,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg,
                                       const std::vector<ModuleBase::Vector3<double>>& kvec_c)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick(" Local_Orbital_wfc","wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int  myid;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info;
    
    //calculate maxnloc for bcasting 2d-wfc
    long maxnloc; // maximum number of elements in local matrix
    info=MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info=MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix

    std::complex<double> **ctot;
    if (out_wfc_lcao && myid == 0)
    {
        ctot = new std::complex<double>*[GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
        }
        ModuleBase::Memory::record("LOWF::ctot", sizeof(std::complex<double>) * GlobalV::NBANDS * GlobalV::NLOCAL);
    }
    
    int naroc[2]; // maximum number of row or column
    for(int iprow=0; iprow<pv->dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            info=MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if(myid==src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work, inc);
                naroc[0]=pv->nrow;
                naroc[1]=pv->ncol_bands;
            }
            info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);
            
            if (out_wfc_lcao)
                info = this->set_wfc_grid(naroc, pv->nb,
                    pv->dim0, pv->dim1, iprow, ipcol,
                    work, wfc_grid, myid, ctot);
            else
                // mohan update 2021-02-12, delte BFIELD option
                info = this->set_wfc_grid(naroc, pv->nb,
                        pv->dim0, pv->dim1, iprow, ipcol,
                        work, wfc_grid);
        }//loop ipcol
    }//loop iprow

    if (out_wfc_lcao && myid == 0)
    {
        std::stringstream ss;
        if (GlobalV::out_app_flag)
        {
            ss << GlobalV::global_out_dir << "LOWF_K_" << ik + 1 << ".dat";
        }
        else
        {
            ss << GlobalV::global_out_dir << istep << "_"
               << "LOWF_K_" << ik + 1 << ".dat";
        }
        ModuleIO::write_wfc_nao_complex(ss.str(), ctot, ik, kvec_c[ik], ekb, wg);
        for (int i = 0; i < GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }
    delete[] work;
    ModuleBase::timer::tick(" Local_Orbital_wfc","wfc_2d_to_grid");
}
#endif