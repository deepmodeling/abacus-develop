#include "grid_technique.h"

#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hsolver/kernels/cuda/helper_cuda.h"
Grid_Technique::Grid_Technique()
{
    this->nlocdimg = nullptr;	
	this->nlocstartg = nullptr;
	this->nad = nullptr;
    this->how_many_atoms = nullptr;
	this->start_ind = nullptr;
	this->which_atom = nullptr;
	this->which_bigcell = nullptr;
	this->which_unitcell = nullptr;
	this->bcell_start = nullptr;
	this->in_this_processor = nullptr;
	this->trace_lo = nullptr;

	this->total_atoms_on_grid = 0;
    allocate_find_R2 = false;
}


Grid_Technique::~Grid_Technique()
{
    delete[] nlocdimg;
    delete[] nlocstartg;
    delete[] nad;
    delete[] how_many_atoms;
	delete[] start_ind;
	delete[] which_atom;
	delete[] which_bigcell;
	delete[] which_unitcell;
	delete[] bcell_start;
	delete[] in_this_processor;
	delete[] trace_lo;
    
    if (allocate_find_R2)
	{
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2[iat];
			delete[] find_R2st[iat];
			delete[] find_R2_sorted_index[iat];
		}
		delete[] find_R2;
		delete[] find_R2st;
		delete[] find_R2_sorted_index;
	}

#if ((defined __CUDA) /* || (defined __ROCM) */)
    for (int i = 0; i < nstreams; ++i)
        checkCudaErrors(cudaStreamDestroy(streams[i]));

    checkCudaErrors(cudaFree(ucell_atom_nwl_g));
    checkCudaErrors(cudaFree(psi_u_g));
    checkCudaErrors(cudaFree(atom_iw2_new_g));
    checkCudaErrors(cudaFree(atom_iw2_ylm_g));
    checkCudaErrors(cudaFree(atom_nw_g));

    checkCudaErrors(cudaFreeHost(psi_input_double_global));
    checkCudaErrors(cudaFreeHost(psi_input_int_global));
    checkCudaErrors(cudaFreeHost(num_psir_global));

    checkCudaErrors(cudaFree(psi_input_double_global_g));
    checkCudaErrors(cudaFree(psi_input_int_global_g));
    checkCudaErrors(cudaFree(num_psir_global_g));
    checkCudaErrors(cudaFree(psir_ylm_left_global_g));
    checkCudaErrors(cudaFree(psir_ylm_right_global_g));

    checkCudaErrors(cudaFreeHost(atom_pair_left_info_global));
    checkCudaErrors(cudaFree(atom_pair_left_info_global_g));

    checkCudaErrors(cudaFreeHost(atom_pair_right_info_global));
    checkCudaErrors(cudaFree(atom_pair_right_info_global_g));

    checkCudaErrors(cudaFreeHost(atom_pair_lda_global));
    checkCudaErrors(cudaFree(atom_pair_lda_global_g));

    checkCudaErrors(cudaFreeHost(atom_pair_ldb_global));
    checkCudaErrors(cudaFree(atom_pair_ldb_global_g));

    checkCudaErrors(cudaFreeHost(atom_pair_ldc_global));
    checkCudaErrors(cudaFree(atom_pair_ldc_global_g));


    checkCudaErrors(cudaFreeHost(atom_pair_left_global));
    checkCudaErrors(cudaFreeHost(atom_pair_right_global));
    checkCudaErrors(cudaFreeHost(atom_pair_output_global));

    checkCudaErrors(cudaFree(atom_pair_left_global_g));
    checkCudaErrors(cudaFree(atom_pair_right_global_g));
    checkCudaErrors(cudaFree(atom_pair_output_global_g));

	const int max_atom_pair_number = GlobalC::ucell.nat * GlobalC::ucell.nat;
    for (int i = 0; i < max_atom_pair_number; i++)
    {
        if (GridVlocal_v2_g[i] != nullptr)
        {
            checkCudaErrors(cudaFree(GridVlocal_v2_g[i]));
        }
    }
    checkCudaErrors(cudaFreeHost(GridVlocal_v2_g));

#endif

}


// This function is called in esolver_ks_lcao_elec.cpp
// after the orbital information has been read,
// this function control the routinue to generate 
// grid technique parameters.
void Grid_Technique::set_pbc_grid(
		const int &ncx_in,
		const int &ncy_in,
		const int &ncz_in,
		const int &bx_in,
		const int &by_in,
		const int &bz_in,
		const int &nbx_in,
		const int &nby_in,
		const int &nbz_in,
		const int &nbxx_in,
		const int &nbzp_start_in,
		const int &nbzp_in,
        const int& ny,
        const int& nplane,
        const int& startz_current)
{
	ModuleBase::TITLE("Grid_Technique","init");
	ModuleBase::timer::tick("Grid_Technique","init");

	if(GlobalV::OUT_LEVEL != "m") 
	{
		GlobalV::ofs_running << "\n SETUP EXTENDED REAL SPACE GRID FOR GRID INTEGRATION" << std::endl;
	}

	// (1) init_meshcell cell and big cell.
	this->set_grid_dim(
		ncx_in,ncy_in,ncz_in,
		bx_in,by_in,bz_in,
		nbx_in,nby_in,nbz_in,
		nbxx_in,nbzp_start_in,nbzp_in);

	this->init_latvec();

	this->init_big_latvec();

	this->init_meshcell_pos();

	// (2) expand the grid
	this->init_grid_expansion();

	// (3) calculate the extended grid.
	this->cal_extended_cell(this->dxe, this->dye, this->dze);

	this->init_tau_in_bigcell();

	// init meshball
	this->delete_meshball_positions(); //LiuXh add 2018-12-14

	this->init_meshball();

	this->init_atoms_on_grid(ny, nplane, startz_current);	

	this->cal_trace_lo();
#if ((defined __CUDA) /* || (defined __ROCM) */)

	this->init_gpu_gint_variables();
#endif

	ModuleBase::timer::tick("Grid_Technique","init");
	return;
}

void Grid_Technique::get_startind(const int& ny, const int& nplane, const int& startz_current)
{
	ModuleBase::TITLE("Grid_Technique","get_startind");

	assert(nbxx>=0);
	// calculates start_ind, which stores the 
	// starting index of each bigcell

	delete[] this->start_ind;
	if(nbxx > 0)
	{
		this->start_ind = new int[nbxx];
		ModuleBase::Memory::record("GT::start_ind", sizeof(int) * nbxx);
		ModuleBase::GlobalFunc::ZEROS(start_ind, nbxx);	
	}
	else
	{
		this->start_ind = nullptr;
		return;
	}

	for(int i=0;i<nbxx;i++)
	{
		int ibx, iby, ibz;
		int ix, iy, iz;

		ibx = i / ( nby * nbzp );
		iby = ( i - ibx * nby * nbzp ) / nbzp;
		ibz = i % nbzp;
	
		ix = ibx * this->bx;
		iy = iby * this->by;
		iz = (ibz + nbzp_start) * this->bz - startz_current;

		int ind = iz + iy * nplane + ix * ny*nplane;
		
		start_ind[i] = ind;
	}

	return;
}

// PLEASE update this 'init_atoms_on_grid' to make
// it adapted to 'cuboid' shape of grid
// mohan add 2021-04-06
void Grid_Technique::init_atoms_on_grid(const int& ny, const int& nplane, const int& startz_current)
{
	ModuleBase::TITLE("Grid_Technique","init_atoms_on_grid");

	assert(nbxx>=0);
	this->get_startind(ny, nplane, startz_current);
	
	// (1) prepare data. 
	// counting the number of atoms whose orbitals have
	// values on the bigcell.
	delete[] this->how_many_atoms;
	if(nbxx > 0)
	{
		this->how_many_atoms = new int[nbxx];
		ModuleBase::Memory::record("GT::how_many_atoms", sizeof(int) * nbxx);
		ModuleBase::GlobalFunc::ZEROS(how_many_atoms, nbxx);
	}	
	else
	{
		this->how_many_atoms = nullptr;
	}
	
	// (2) information about gloabl grid
	// and local grid.
	// mohan add 2010-07-02
	int *ind_bigcell;
	bool *bigcell_on_processor; // normal local form.
	this->check_bigcell(ind_bigcell, bigcell_on_processor);

	// (3) Find the atoms using
	// when doing grid integration. 
	delete[] in_this_processor;
	this->in_this_processor = new bool[GlobalC::ucell.nat];
	for(int i=0; i<GlobalC::ucell.nat; i++)
	{
		in_this_processor[i] = false;
	}

	// init atoms on grid	
	assert( this->nxyze > 0);
	int* index2normal = new int[this->nxyze];
	assert( index2normal != NULL );
	ModuleBase::Memory::record("GT::index2normal", sizeof(int) * this->nxyze);
	this->grid_expansion_index(1,index2normal); 

	// (5) record how many atoms on
	// each local grid point (ix,iy,iz)
	int iat=0;
	int normal;
	this->total_atoms_on_grid = 0;
	int nat_local = 0;
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			for(int im=0; im< this->meshball_ncells; im++)
			{
				// bcell[iat]: which bcell iat atom is in.
				// ball[im]: relative position of adjacent bcell.
				normal = index2normal[ this->index_atom[iat] + this->index_ball[im] ];

				if(normal >= nbxyz)
				{
					std::cout << " index_atom=" << index_atom[iat] << std::endl;
					std::cout << " index_ball=" << index_ball[im] << std::endl;
					std::cout << " normal=" << normal << std::endl;
					std::cout << " nbxyz=" << nbxyz << std::endl;
					ModuleBase::WARNING_QUIT("Grid_Technique::init_atoms_on_grid","normal >= nbxyz");
				}

				assert(normal>=0);

				int f = ind_bigcell[normal];
				if(!bigcell_on_processor[normal]) continue;

				++how_many_atoms[f];
				++total_atoms_on_grid;

				this->in_this_processor[iat] = true;
			}
			if(this->in_this_processor[iat]) ++nat_local;
			++iat;
		}
	}

	delete[] ind_bigcell;
	delete[] bigcell_on_processor;

	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total_atoms_on_grid",total_atoms_on_grid);
	
	int stop = 0;
	if(total_atoms_on_grid == 0)
	{
		GlobalV::ofs_running << " No atoms on this sub-FFT-mesh." << std::endl;
		stop = 1;
	} 	
    Parallel_Reduce::reduce_all(stop);
	if(stop)
	{
		ModuleBase::WARNING("Grid_Technique::init_atoms_on_grid","No atom on this sub-FFT-mesh.");
	}

	// calculate the trach of local ia to global iat
	if(nat_local>0)
	{
		this->trace_iat.resize(nat_local);
		for(int iat=GlobalC::ucell.nat-1; iat>=0;iat--)
		{
			if(this->in_this_processor[iat])
			{
				this->trace_iat[--nat_local] = iat;
			}
		}
	}

	// need how_many_atoms first.
	this->cal_grid_integration_index();
	// bcell_start is needed.
	this->init_atoms_on_grid2(index2normal);
	delete[] index2normal;	
	return;
}

void Grid_Technique::check_bigcell(int* &ind_bigcell, bool* &bigcell_on_processor)
{
	//check if a given bigcell is treated on this processor
	const int zstart = nbzp_start;
	const int zend = nbzp + zstart;
	const int nbyz = nby * nbz;
	const int nz = nbzp;

	int iz_now, ix, iy, iz, ind;
	bool flag;

	ind_bigcell = new int[nbxyz];
	bigcell_on_processor=new bool[nbxyz];
	for(int i=0;i<nbxyz;i++)
	{
		int iz_now = i % nbz;
		if(iz_now<zstart || iz_now>=zend)
		{
			flag=false;
		}
		else
		{
			flag=true;
			ix = i / nbyz;
			iy = ( i - ix * nbyz ) / nbz;
			iz = iz_now - zstart;
			ind = ix * nby * nz + iy * nz + iz;
			//no need to calculate index if bigcell is
			//not on this processor
		}

		ind_bigcell[i]=ind;
		bigcell_on_processor[i]=flag;
	}			
	return;
}

void Grid_Technique::init_atoms_on_grid2(const int* index2normal)
{	
	ModuleBase::TITLE("Grid_Techinique","init_atoms_on_grid2");

	if(total_atoms_on_grid==0) 
	{	
		ModuleBase::WARNING("Grid_Technique::init_atoms_on_grid2","no atom on this sub FFT grid.");
		return;
	}

	int* index2ucell = new int[this->nxyze];
	assert( index2ucell != NULL );
	ModuleBase::Memory::record("GT::index2ucell", sizeof(int) * this->nxyze);	
	this->grid_expansion_index(0,index2ucell);
	
	int *ind_bigcell;
	bool *bigcell_on_processor; // normal local form.
	this->check_bigcell(ind_bigcell, bigcell_on_processor);

	//--------------------------------------
	// save which atom is in the bigcell.
	//--------------------------------------
	delete[] which_atom;
	this->which_atom = new int[total_atoms_on_grid];
	assert( which_atom != 0);
	ModuleBase::Memory::record("GT::which_atom", sizeof(int) * total_atoms_on_grid);

	delete[] which_bigcell;
	this->which_bigcell = new int[total_atoms_on_grid];
	assert( which_bigcell != 0);
	ModuleBase::Memory::record("GT::which_bigcell", sizeof(int) * total_atoms_on_grid);

	delete[] which_unitcell;
	this->which_unitcell = new int[total_atoms_on_grid];
	assert( which_unitcell != 0);
	ModuleBase::Memory::record("GT::which_unitcell", sizeof(int) * total_atoms_on_grid);
	// for each atom, first we need to locate which cell
	// the atom is in, then we search meshball aroung this
	// grid, and record each grid's atom position.
	int count = 0;
	int iat = 0;
	ModuleBase::GlobalFunc::ZEROS(this->how_many_atoms, nbxx);
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			// zero bigcell of meshball indicate ?
			for(int im=0; im< this->meshball_ncells; im++)
			{
				const int extgrid = this->index_atom[iat] + this->index_ball[im];
				const int normal = index2normal[ extgrid ];
			
				// mohan add 2010-07-01
				int f = ind_bigcell[normal];
				if(!bigcell_on_processor[normal]) continue;
				
				// it's not the normal order to calculate which_atom
				// and which_bigcell, especailly in 1D array. 
				// Each grid's adjacent atom number is different,
				// so, first we need to locate which grid, using
				// bcell_start, then we need to count which adjacent atom.
				// using how_many_atoms.
				int index = this->bcell_start[f] + this->how_many_atoms[f];
				
				// we save which_atom and which_bigcell in 1D array,
				// once you want to use this in grid integration, 
				// the only information you got is the 'normal' index, 
				// so you need to use bcell_start
				// to get the 'mesh_index', then you can you this mesh_index
				// to use which_atom or which_bigcell.
				this->which_atom[ index ] = iat;
				this->which_bigcell[ index ] = im;
				this->which_unitcell[ index ] = index2ucell[extgrid];
				 
				++how_many_atoms[f];
				++count;
			}
			++iat;
		}
	}
	assert( count == total_atoms_on_grid );
	delete[] index2ucell;
	delete[] ind_bigcell;
	delete[] bigcell_on_processor;
	return;
}

void Grid_Technique::cal_grid_integration_index(void)
{
	// save the start 
	delete[] this->bcell_start;
	if(nbxx > 0)
	{
		this->bcell_start = new int[nbxx];
		ModuleBase::Memory::record("GT::bcell_start", sizeof(int) * nbxx);
		this->bcell_start[0] = 0;
		for(int i=1; i<nbxx; i++)
		{
			this->bcell_start[i] = this->bcell_start[i-1] + this->how_many_atoms[i-1];
		}
	}
	else
	{
		this->bcell_start = nullptr;
	}
	// calculate which grid has the largest number of atoms,
	// and how many atoms.
	this->max_atom = 0;
	for(int i=0; i<nbxx; i++)
	{
		this->max_atom = std::max( this->max_atom, this->how_many_atoms[i]);
	}

#ifdef __MPI
	int* all = new int[GlobalV::NPROC];
	ModuleBase::GlobalFunc::ZEROS(all, GlobalV::NPROC);
	Parallel_Reduce::gather_int_all(max_atom,all);
	if(GlobalV::MY_RANK==0)
	{
		GlobalV::ofs_warning << std::setw(15) << "Processor" << std::setw(15) << "Atom" << std::endl;
		for(int i=0; i<GlobalV::NPROC; i++)
		{
			GlobalV::ofs_warning << std::setw(15) << i+1 << std::setw(15) << all[i] << std::endl; 
		}
	}
	delete[] all;
#endif

	if(GlobalV::test_gridt)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Max atom on bigcell",max_atom);
	return;
}

// set 'lgd' variable
void Grid_Technique::cal_trace_lo(void)
{	
	ModuleBase::TITLE("Grid_Technique","cal_trace_lo");
	// save the atom information in trace_lo,
	// in fact the trace_lo dimension can be reduced
	// to GlobalC::ucell.nat, but I think this is another way.
	delete[] trace_lo;
	this->trace_lo = new int[GlobalV::NLOCAL];
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		this->trace_lo[i] = -1;
	}
	ModuleBase::Memory::record("GT::trace_lo", sizeof(int) * GlobalV::NLOCAL);

	this->lnat = 0;
	this->lgd = 0;
	int iat = 0;
	int iw_all=0;
	int iw_local=0;

	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			if(this->in_this_processor[iat])
			{
				++lnat;
				int nw0 = GlobalC::ucell.atoms[it].nw;
				if(GlobalV::NSPIN==4)
				{//added by zhengdy-soc, need to be double in soc
					nw0 *= 2;
					this->lgd += nw0;
				}
				else
				{
					this->lgd += GlobalC::ucell.atoms[it].nw;
				}
				
				for(int iw=0; iw<nw0; iw++)
				{
					this->trace_lo[iw_all] = iw_local;
					++iw_local; 
					++iw_all;
				}
			}
			else
			{
				// global index of atomic orbitals
				iw_all += GlobalC::ucell.atoms[it].nw;
				if(GlobalV::NSPIN==4) iw_all += GlobalC::ucell.atoms[it].nw;
			}
			++iat;
		}
	}
	
	//------------
	// for test
	//------------
//	for(int i=0; i<GlobalV::NLOCAL; ++i)
//	{
//		GlobalV::ofs_running << " i=" << i+1 << " trace_lo=" << trace_lo[i] << std::endl;
//	}

	if(GlobalV::OUT_LEVEL != "m") 
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Atom number in sub-FFT-grid",lnat);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Local orbitals number in sub-FFT-grid",lgd);
	}

	assert(iw_local == lgd);
	assert(iw_all == GlobalV::NLOCAL);
	return;
}


#if ((defined __CUDA) /* || (defined __ROCM) */)

	void Grid_Technique::init_gpu_gint_variables()
	{
		double ylmcoef[100];
		ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
		for (int i = 0; i < 100; i++)
		{
			ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
		}
		checkCudaErrors(cudaMalloc((void **)&ylmcoef_g, 100 * sizeof(double)));
		checkCudaErrors(cudaMemcpy(ylmcoef_g, ylmcoef, 100 * sizeof(double), cudaMemcpyHostToDevice));

		const Numerical_Orbital_Lm *pointer;
		double max_cut = 0;
		for (int i = 0; i < GlobalC::ucell.ntype; i++)
		{
			if (GlobalC::ORB.Phi[i].getRcut() > max_cut)
			{
				max_cut = GlobalC::ORB.Phi[i].getRcut();
			}
		}

		int atom_nw_now[GlobalC::ucell.ntype];
		int ucell_atom_nwl_now[GlobalC::ucell.ntype];
		for (int i = 0; i < GlobalC::ucell.ntype; i++)
		{
			atom_nw_now[i] = GlobalC::ucell.atoms[i].nw;
			ucell_atom_nwl_now[i] = GlobalC::ucell.atoms[i].nwl;
		}

		nr_max = static_cast<int>(1000 * max_cut) + 10;
		double psi_u_now[GlobalC::ucell.ntype * GlobalC::ucell.nwmax * nr_max * 2];
		memset(psi_u_now, 0, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * nr_max * 2 * sizeof(double));
		bool atom_iw2_new_now[GlobalC::ucell.ntype * GlobalC::ucell.nwmax];
		memset(atom_iw2_new_now, 0, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(bool));
		int atom_iw2_ylm_now[GlobalC::ucell.ntype * GlobalC::ucell.nwmax];
		memset(atom_iw2_ylm_now, 0, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(int));

		Atom *atomx;
		for (int i = 0; i < GlobalC::ucell.ntype; i++)
		{
			atomx = &GlobalC::ucell.atoms[i];
			for (int j = 0; j < GlobalC::ucell.nwmax; j++)
			{
				if (j < atomx->nw)
				{
					atom_iw2_new_now[i * GlobalC::ucell.nwmax + j] = atomx->iw2_new[j];
					atom_iw2_ylm_now[i * GlobalC::ucell.nwmax + j] = atomx->iw2_ylm[j];
					pointer = &GlobalC::ORB.Phi[i].PhiLN(atomx->iw2l[j], atomx->iw2n[j]);
					for (int k = 0; k < nr_max; k++)
					{
						int index_temp = (i * GlobalC::ucell.nwmax * nr_max + j * nr_max + k) * 2;
						if (k < pointer->nr_uniform)
						{
							psi_u_now[index_temp] = pointer->psi_uniform[k];
							psi_u_now[index_temp + 1] = pointer->dpsi_uniform[k];
						}
					}
				}
			}
		}

		checkCudaErrors(cudaMalloc((void **)&atom_nw_g, GlobalC::ucell.ntype * sizeof(int)));
		checkCudaErrors(cudaMemcpy(atom_nw_g, atom_nw_now, GlobalC::ucell.ntype * sizeof(int), cudaMemcpyHostToDevice));

		checkCudaErrors(cudaMalloc((void **)&ucell_atom_nwl_g, GlobalC::ucell.ntype * sizeof(int)));
		checkCudaErrors(cudaMemcpy(ucell_atom_nwl_g, ucell_atom_nwl_now, GlobalC::ucell.ntype * sizeof(int), cudaMemcpyHostToDevice));

		checkCudaErrors(cudaMalloc((void **)&psi_u_g, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * nr_max * sizeof(double) * 2));
		checkCudaErrors(cudaMemcpy(psi_u_g, psi_u_now, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * nr_max * sizeof(double) * 2, cudaMemcpyHostToDevice));

		checkCudaErrors(cudaMalloc((void **)&atom_iw2_new_g, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(bool)));
		checkCudaErrors(cudaMalloc((void **)&atom_iw2_ylm_g, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(int)));
		checkCudaErrors(cudaMemcpy(atom_iw2_new_g, atom_iw2_new_now, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(bool), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(atom_iw2_ylm_g, atom_iw2_ylm_now, GlobalC::ucell.ntype * GlobalC::ucell.nwmax * sizeof(int), cudaMemcpyHostToDevice));

		const int max_atom_pair_number = GlobalC::ucell.nat * GlobalC::ucell.nat;
		checkCudaErrors(cudaMallocHost((void **)&GridVlocal_v2_g, max_atom_pair_number * sizeof(double *)));  // the points to gpu memory, but gpu memory address save on host

		for (int iat1 = 0; iat1 < GlobalC::ucell.nat; iat1++)
		{
			int it1 = GlobalC::ucell.iat2it[iat1];
			int lo1 = trace_lo[GlobalC::ucell.itiaiw2iwt(it1, GlobalC::ucell.iat2ia[iat1], 0)];
			for (int iat2 = 0; iat2 < GlobalC::ucell.nat; iat2++)
			{
				int it2 = GlobalC::ucell.iat2it[iat2];
				int lo2 = trace_lo[GlobalC::ucell.itiaiw2iwt(it2, GlobalC::ucell.iat2ia[iat2], 0)];
				if (lo1 <= lo2) {
					int atom_pair_nw = GlobalC::ucell.atoms[it1].nw * GlobalC::ucell.atoms[it2].nw;
					checkCudaErrors(cudaMalloc((void **)&GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2], atom_pair_nw * sizeof(double)));
					checkCudaErrors(cudaMemset(GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2], 0, atom_pair_nw * sizeof(double)));
				} else {
					GridVlocal_v2_g[iat1 * GlobalC::ucell.nat + iat2] = nullptr;
				}
			}
		}

		psir_size = nbzp * max_atom * bxyz * GlobalC::ucell.nwmax;
		checkCudaErrors(cudaMalloc((void **)&psir_ylm_left_global_g, psir_size * nstreams * sizeof(double)));
		checkCudaErrors(cudaMalloc((void **)&psir_ylm_right_global_g, psir_size * nstreams * sizeof(double)));
		checkCudaErrors(cudaMemset(psir_ylm_left_global_g, 0, psir_size * nstreams * sizeof(double)));
		checkCudaErrors(cudaMemset(psir_ylm_right_global_g, 0, psir_size * nstreams * sizeof(double)));

		atom_pair_size_of_meshcell = max_atom * max_atom;
		atom_pair_size_over_nbz = atom_pair_size_of_meshcell * nbzp;

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_left_info_global, atom_pair_size_over_nbz * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_left_info_global_g, atom_pair_size_over_nbz * nstreams * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_right_info_global, atom_pair_size_over_nbz * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_right_info_global_g, atom_pair_size_over_nbz * nstreams * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_lda_global, atom_pair_size_over_nbz * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_lda_global_g, atom_pair_size_over_nbz * nstreams * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_ldb_global, atom_pair_size_over_nbz * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_ldb_global_g, atom_pair_size_over_nbz * nstreams * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_ldc_global, atom_pair_size_over_nbz * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_ldc_global_g, atom_pair_size_over_nbz * nstreams * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&atom_pair_left_global, atom_pair_size_over_nbz * nstreams * sizeof(double *)));
		checkCudaErrors(cudaMallocHost((void **)&atom_pair_right_global, atom_pair_size_over_nbz * nstreams * sizeof(double *)));
		checkCudaErrors(cudaMallocHost((void **)&atom_pair_output_global, atom_pair_size_over_nbz * nstreams * sizeof(double *)));

		checkCudaErrors(cudaMalloc((void **)&atom_pair_left_global_g, atom_pair_size_over_nbz * nstreams * sizeof(double *)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_right_global_g, atom_pair_size_over_nbz * nstreams * sizeof(double *)));
		checkCudaErrors(cudaMalloc((void **)&atom_pair_output_global_g, atom_pair_size_over_nbz * nstreams * sizeof(double *)));

		psi_size_max = max_atom * bxyz * nbzp;
		psi_size_max_per_z = max_atom * bxyz;
		checkCudaErrors(cudaMallocHost((void **)&psi_input_double_global, psi_size_max * nstreams * 5 * sizeof(double)));
		checkCudaErrors(cudaMalloc((void **)&psi_input_double_global_g, psi_size_max * nstreams * 5 * sizeof(double)));

		checkCudaErrors(cudaMallocHost((void **)&psi_input_int_global, psi_size_max * nstreams * 2 * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&psi_input_int_global_g, psi_size_max * nstreams * 2 * sizeof(int)));

		checkCudaErrors(cudaMallocHost((void **)&num_psir_global, nbzp * nstreams * sizeof(int)));
		checkCudaErrors(cudaMalloc((void **)&num_psir_global_g, nbzp * nstreams * sizeof(int)));

		for (int i = 0; i < nstreams; ++i)
		{
			checkCudaErrors(cudaStreamCreate(&streams[i]));
		}

		gemm_algo_selector(GlobalC::ucell.nwmax,GlobalC::ucell.nwmax, bxyz, fastest_matrix_mul);
	}
#endif
