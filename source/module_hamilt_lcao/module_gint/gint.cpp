#include "gint.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "gint_gpu/gint_vl.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint::~Gint()
{
	delete this->hRGint;
	delete this->hRGintCd;
	
}

void Gint::cal_gint(Gint_inout *inout)
{

	ModuleBase::timer::tick("Gint_interface", "cal_gint");

	if(inout->job==Gint_Tools::job_type::vlocal) ModuleBase::TITLE("Gint_interface","cal_gint_vlocal");
	if(inout->job==Gint_Tools::job_type::vlocal_meta) ModuleBase::TITLE("Gint_interface","cal_gint_vlocal_meta");
	if(inout->job==Gint_Tools::job_type::rho) ModuleBase::TITLE("Gint_interface","cal_gint_rho");
	if(inout->job==Gint_Tools::job_type::tau) ModuleBase::TITLE("Gint_interface","cal_gint_tau");
	if(inout->job==Gint_Tools::job_type::force) ModuleBase::TITLE("Gint_interface","cal_gint_force");
	if(inout->job==Gint_Tools::job_type::force_meta) ModuleBase::TITLE("Gint_interface","cal_gint_force_meta");

	if(inout->job==Gint_Tools::job_type::vlocal) ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
	if(inout->job==Gint_Tools::job_type::vlocal_meta) ModuleBase::timer::tick("Gint_interface","cal_gint_vlocal_meta");
	if(inout->job==Gint_Tools::job_type::rho) ModuleBase::timer::tick("Gint_interface","cal_gint_rho");
	if(inout->job==Gint_Tools::job_type::tau) ModuleBase::timer::tick("Gint_interface","cal_gint_tau");
	if(inout->job==Gint_Tools::job_type::force) ModuleBase::timer::tick("Gint_interface","cal_gint_force");
	if(inout->job==Gint_Tools::job_type::force_meta) ModuleBase::timer::tick("Gint_interface","cal_gint_force_meta");

	const int max_size = this->gridt->max_atom;
	const int LD_pool = max_size*GlobalC::ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;

    if(max_size!=0)
    {
		if (inout->job == Gint_Tools::job_type::vlocal && GlobalV::GAMMA_ONLY_LOCAL && lgd > 0)
		{
			double ylmcoef[100];
			ModuleBase::GlobalFunc::ZEROS(ylmcoef, 100);
			for (int i = 0; i < 100; i++)
			{
				ylmcoef[i] = ModuleBase::Ylm::ylmcoef[i];
			}
			gint_gamma_vl_gpu(this->hRGint,
							lgd,
							max_size,
							GlobalC::ucell.omega / this->ncxyz,
							inout->vl,
							ylmcoef,
							this->nplane,
							GlobalV::NLOCAL,
							this->nbxx,
							*this->gridt);
			ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
			return;
		}
		else
		{
#ifdef __MKL
		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(1);
#endif

#ifdef _OPENMP
    	#pragma omp parallel
#endif
		{
            //prepare some constants
			const int ncyz = this->ny*this->nplane; // mohan add 2012-03-25
			const double dv = GlobalC::ucell.omega/this->ncxyz;

			// it's a uniform grid to save orbital values, so the delta_r is a constant.
			const double delta_r = GlobalC::ORB.dr_uniform;

            if((inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta) && !GlobalV::GAMMA_ONLY_LOCAL)
            {
                if(!pvpR_alloc_flag)
                {
                    ModuleBase::WARNING_QUIT("Gint_interface::cal_gint","pvpR has not been allocated yet!");
                }
                else
                {
                    ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin], nnrg);
                }
            }

			if(inout->job==Gint_Tools::job_type::dvlocal)
			{
				if(GlobalV::GAMMA_ONLY_LOCAL)
				{
					ModuleBase::WARNING_QUIT("Gint_interface::cal_gint","dvlocal only for k point!");
				}
				ModuleBase::GlobalFunc::ZEROS(this->pvdpRx_reduced[inout->ispin], nnrg);
				ModuleBase::GlobalFunc::ZEROS(this->pvdpRy_reduced[inout->ispin], nnrg);
				ModuleBase::GlobalFunc::ZEROS(this->pvdpRz_reduced[inout->ispin], nnrg);
			}

            //perpare auxiliary arrays to store thread-specific values
#ifdef _OPENMP
			double* pvpR_thread;
			hamilt::HContainer<double>* hRGint_thread = nullptr;// auxiliary pointer for multi-threading
			if(inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta)
			{
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    pvpR_thread = new double[nnrg];
                    ModuleBase::GlobalFunc::ZEROS(pvpR_thread, nnrg);
                }
                if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
                {
					hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);
                }
			}

			double *pvdpRx_thread, *pvdpRy_thread, *pvdpRz_thread;
			if(inout->job==Gint_Tools::job_type::dvlocal)
			{
				pvdpRx_thread = new double[nnrg];
				ModuleBase::GlobalFunc::ZEROS(pvdpRx_thread, nnrg);
				pvdpRy_thread = new double[nnrg];
				ModuleBase::GlobalFunc::ZEROS(pvdpRy_thread, nnrg);
				pvdpRz_thread = new double[nnrg];
				ModuleBase::GlobalFunc::ZEROS(pvdpRz_thread, nnrg);								
			}

			ModuleBase::matrix fvl_dphi_thread;
			ModuleBase::matrix svl_dphi_thread;
			if(inout->job==Gint_Tools::job_type::force || inout->job==Gint_Tools::job_type::force_meta)
			{
				if(inout->isforce)
				{
					fvl_dphi_thread.create(inout->fvl_dphi->nr,inout->fvl_dphi->nc);
					fvl_dphi_thread.zero_out();
				}
				if(inout->isstress)
				{
					svl_dphi_thread.create(inout->svl_dphi->nr,inout->svl_dphi->nc);
					svl_dphi_thread.zero_out();
				}
			}

    		#pragma omp for
#endif
            // entering the main loop of grid points
			for(int grid_index = 0; grid_index < this->nbxx; grid_index++)
			{
				// get the value: how many atoms has orbital value on this grid.
				const int na_grid = this->gridt->how_many_atoms[ grid_index ];

				if(na_grid==0) continue;

				if(inout->job == Gint_Tools::job_type::rho)
				{
					//int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
                    int* vindex = Gint_Tools::get_vindex(this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz);
                    this->gint_kernel_rho(na_grid, grid_index, delta_r, vindex, LD_pool, inout);
					delete[] vindex;
				}
				else if(inout->job == Gint_Tools::job_type::tau)
				{
                    int* vindex = Gint_Tools::get_vindex(this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz);
                    this->gint_kernel_tau(na_grid, grid_index, delta_r, vindex, LD_pool, inout);
					delete[] vindex;
				}
				else if(inout->job == Gint_Tools::job_type::force)
				{
                    double* vldr3 = Gint_Tools::get_vldr3(inout->vl, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
                    double** DM_in;
					if(GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM[GlobalV::CURRENT_SPIN];
					if(!GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM_R;
					#ifdef _OPENMP
						this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool,
							DM_in, inout->ispin, inout->isforce, inout->isstress,
							&fvl_dphi_thread, &svl_dphi_thread);
					#else
						this->gint_kernel_force(na_grid, grid_index, delta_r, vldr3, LD_pool,
							DM_in, inout->ispin, inout->isforce, inout->isstress,
							inout->fvl_dphi, inout->svl_dphi);
					#endif
					delete[] vldr3;
				}
				else if(inout->job==Gint_Tools::job_type::vlocal)
				{
                    double* vldr3 = Gint_Tools::get_vldr3(inout->vl, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
#ifdef _OPENMP
						if((GlobalV::GAMMA_ONLY_LOCAL && lgd>0) || !GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
								pvpR_thread, hRGint_thread);
						}
					#else
						if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool, pvpR_grid);
						}
						if(!GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
								this->pvpR_reduced[inout->ispin]);
						}
					#endif
					delete[] vldr3;
				}
				else if(inout->job==Gint_Tools::job_type::dvlocal)
				{
                    double* vldr3 = Gint_Tools::get_vldr3(inout->vl, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
#ifdef _OPENMP
						this->gint_kernel_dvlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
							pvdpRx_thread, pvdpRy_thread, pvdpRz_thread);
					#else
						this->gint_kernel_dvlocal(na_grid, grid_index, delta_r, vldr3, LD_pool,
							this->pvdpRx_reduced[inout->ispin], this->pvdpRy_reduced[inout->ispin], this->pvdpRz_reduced[inout->ispin]);
					#endif
					delete[] vldr3;
				}
				else if(inout->job==Gint_Tools::job_type::vlocal_meta)
				{
                    double* vldr3 = Gint_Tools::get_vldr3(inout->vl, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
                    double* vkdr3 = Gint_Tools::get_vldr3(inout->vofk,this->bxyz,  this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
#ifdef _OPENMP
						if((GlobalV::GAMMA_ONLY_LOCAL && lgd>0) || !GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
								pvpR_thread, hRGint_thread);
						}
					#else
						if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool, pvpR_grid);
						}
						if(!GlobalV::GAMMA_ONLY_LOCAL)
						{
							this->gint_kernel_vlocal_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
								this->pvpR_reduced[inout->ispin]);
						}
					#endif
					delete[] vldr3;
					delete[] vkdr3;
				}
				else if(inout->job == Gint_Tools::job_type::force_meta)
				{
                    double* vldr3 = Gint_Tools::get_vldr3(inout->vl, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
                    double* vkdr3 = Gint_Tools::get_vldr3(inout->vofk, this->bxyz, this->bx, this->by, this->bz,
                        this->nplane, this->gridt->start_ind[grid_index], ncyz, dv);
                    double** DM_in;
					if(GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM[GlobalV::CURRENT_SPIN];
					if(!GlobalV::GAMMA_ONLY_LOCAL) DM_in = inout->DM_R;
					#ifdef _OPENMP
						this->gint_kernel_force_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
							DM_in, inout->ispin, inout->isforce, inout->isstress,
							&fvl_dphi_thread, &svl_dphi_thread);
					#else
						this->gint_kernel_force_meta(na_grid, grid_index, delta_r, vldr3, vkdr3, LD_pool,
							DM_in, inout->ispin, inout->isforce, inout->isstress,
							inout->fvl_dphi, inout->svl_dphi);
					#endif
					delete[] vldr3;
					delete[] vkdr3;
				}
			} // int grid_index
#ifdef _OPENMP
			if(inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta)
			{
                if(GlobalV::GAMMA_ONLY_LOCAL && lgd>0)
                {
                    #pragma omp critical(gint_gamma)
                    {
						this->hRGint->add(*hRGint_thread);
					}
                    delete hRGint_thread;
                }
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    #pragma omp critical(gint_k)
                    for(int innrg=0; innrg<nnrg; innrg++)
                    {
                        pvpR_reduced[inout->ispin][innrg] += pvpR_thread[innrg];
                    }
                    delete[] pvpR_thread;
                }
			}

			#pragma omp critical(gint)
			if(inout->job==Gint_Tools::job_type::force  || inout->job==Gint_Tools::job_type::force_meta)
			{
				if(inout->isforce)
				{
					inout->fvl_dphi[0]+=fvl_dphi_thread;
				}
				if(inout->isstress)
				{
					inout->svl_dphi[0]+=svl_dphi_thread;
				}
			}
#endif
		} // end of #pragma omp parallel
		#ifdef __MKL
			mkl_set_num_threads(mkl_threads);
		#endif
		}

    } // end of if (max_size)

	ModuleBase::timer::tick("Gint_interface", "cal_gint");

	if(inout->job==Gint_Tools::job_type::vlocal) ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
	if(inout->job==Gint_Tools::job_type::vlocal_meta) ModuleBase::timer::tick("Gint_interface","cal_gint_vlocal_meta");
	if(inout->job==Gint_Tools::job_type::rho) ModuleBase::timer::tick("Gint_interface","cal_gint_rho");
	if(inout->job==Gint_Tools::job_type::tau) ModuleBase::timer::tick("Gint_interface","cal_gint_tau");
	if(inout->job==Gint_Tools::job_type::force) ModuleBase::timer::tick("Gint_interface","cal_gint_force");
	if(inout->job==Gint_Tools::job_type::force_meta) ModuleBase::timer::tick("Gint_interface","cal_gint_force_meta");
	return;
}

void Gint::prep_grid(
    const Grid_Technique& gt,
    const int& nbx_in,
	const int &nby_in,
	const int &nbz_in,
	const int &nbz_start_in,
    const int& ncxyz_in,
    const int& bx_in,
    const int& by_in,
    const int& bz_in,
    const int& bxyz_in,
    const int& nbxx_in,
    const int& ny_in,
    const int& nplane_in,
    const int& startz_current_in)
{
	ModuleBase::TITLE(GlobalV::ofs_running,"Gint_k","prep_grid");

    this->gridt = &gt;
    this->nbx = nbx_in;
	this->nby = nby_in;
	this->nbz = nbz_in;
	this->ncxyz = ncxyz_in;
    this->nbz_start = nbz_start_in;
    this->bx = bx_in;
    this->by = by_in;
    this->bz = bz_in;
    this->bxyz = bxyz_in;
    this->nbxx = nbxx_in;
    this->ny = ny_in;
    this->nplane = nplane_in;
    this->startz_current = startz_current_in;
    assert(nbx > 0);
	assert(nby>0);
	assert(nbz>=0);
    assert(ncxyz > 0);
    assert(bx > 0);
    assert(by > 0);
    assert(bz > 0);
    assert(bxyz > 0);
    assert(nbxx >= 0);
    assert(ny > 0);
    assert(nplane >= 0);
    assert(startz_current >= 0);

	assert( GlobalC::ucell.omega > 0.0);

	return;
}

void Gint::initialize_pvpR(
	const UnitCell& ucell_in,
	Grid_Driver* gd)
{
	ModuleBase::TITLE("Gint","initialize_pvpR");

	int npol = 1;
	// there is the only resize code of DMRGint
	if(this->DMRGint.size() == 0)
	{
		this->DMRGint.resize(GlobalV::NSPIN);
	}
	if(GlobalV::NSPIN!=4)
	{
		if(this->hRGint != nullptr) delete this->hRGint;
		this->hRGint = new hamilt::HContainer<double>(ucell_in.nat);
	}
	else
	{
		npol = 2;
		if(this->hRGintCd != nullptr) delete this->hRGintCd;
		this->hRGintCd = new hamilt::HContainer<std::complex<double>>(ucell_in.nat);
		for (int is = 0; is < GlobalV::NSPIN; is++)
		{
			if (this->DMRGint[is] != nullptr)
			{
				delete this->DMRGint[is];
			}
			this->DMRGint[is] = new hamilt::HContainer<double>(ucell_in.nat);
		}
#ifdef __MPI
		if(this->DMRGint_full != nullptr) delete this->DMRGint_full;
		this->DMRGint_full = new hamilt::HContainer<double>(ucell_in.nat);
#endif
	}

	// prepare the row_index and col_index for construct AtomPairs, they are same, name as orb_index
	std::vector<int> orb_index(ucell_in.nat + 1);
	orb_index[0] = 0;
	for(int i=1;i<orb_index.size();i++)
	{
		int type = ucell_in.iat2it[i-1];
		orb_index[i] = orb_index[i-1] + ucell_in.atoms[type].nw;
	}
	std::vector<int> orb_index_npol;
	if(npol == 2)
	{
		orb_index_npol.resize(ucell_in.nat + 1);
		orb_index_npol[0] = 0;
		for(int i=1;i<orb_index_npol.size();i++)
		{
			int type = ucell_in.iat2it[i-1];
			orb_index_npol[i] = orb_index_npol[i-1] + ucell_in.atoms[type].nw * npol;
		}
	}

	if(GlobalV::GAMMA_ONLY_LOCAL && GlobalV::NSPIN != 4)
	{
		this->hRGint->fix_gamma();
	}
	for (int T1 = 0; T1 < ucell_in.ntype; ++T1)
	{
		const Atom* atom1 = &(ucell_in.atoms[T1]);
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			auto& tau1 = atom1->tau[I1];

			gd->Find_atom(ucell_in, tau1, T1, I1);

			const int iat1 = ucell_in.itia2iat(T1,I1);

			// for grid integration (on FFT box),
			// we only need to consider <phi_i | phi_j>,

			// whether this atom is in this processor.
			if(this->gridt->in_this_processor[iat1])
			{
				for (int ad = 0; ad < gd->getAdjacentNum()+1; ++ad)
				{
					const int T2 = gd->getType(ad);
					const int I2 = gd->getNatom(ad);
					const int iat2 = ucell_in.itia2iat(T2, I2);
					const Atom* atom2 = &(ucell_in.atoms[T2]); 

					// NOTE: hRGint wil save total number of atom pairs, 
					// if only upper triangle is saved, the lower triangle will be lost in 2D-block parallelization.
					// if the adjacent atom is in this processor.
					if(this->gridt->in_this_processor[iat2])
					{
						ModuleBase::Vector3<double> dtau = gd->getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * ucell_in.lat0;
						double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

						//if(distance < rcut)
						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
							// calculate R index
							auto& R_index = gd->getBox(ad);
							// insert this atom-pair into this->hRGint
							if(npol == 1)
							{
								hamilt::AtomPair<double> tmp_atom_pair(iat1, iat2, R_index.x, R_index.y, R_index.z, orb_index.data(), orb_index.data(), ucell_in.nat);
								this->hRGint->insert_pair(tmp_atom_pair);
							}
							else
							{
								// HR is complex and size is nw * npol
								hamilt::AtomPair<std::complex<double>> tmp_atom_pair(iat1, iat2, R_index.x, R_index.y, R_index.z, orb_index_npol.data(), orb_index_npol.data(), ucell_in.nat);
								this->hRGintCd->insert_pair(tmp_atom_pair);
								// DMR is double now and size is nw
								hamilt::AtomPair<double> tmp_dmR(iat1, iat2, R_index.x, R_index.y, R_index.z, orb_index.data(), orb_index.data(), ucell_in.nat);
								for (int is = 0; is < this->DMRGint.size(); is++)
								{
									this->DMRGint[is]->insert_pair(tmp_dmR);
								}
#ifdef __MPI					
								hamilt::AtomPair<double> tmp_dmR_full(iat1, iat2, R_index.x, R_index.y, R_index.z, orb_index_npol.data(), orb_index_npol.data(), ucell_in.nat);
								// tmp DMR for transfer
								this->DMRGint_full->insert_pair(tmp_dmR_full);
#endif
							}
						}
					}// end iat2
				}// end ad
			}// end iat
		}// end I1
	}// end T1
	if(npol == 1)
	{
		this->hRGint->allocate(0);
		ModuleBase::Memory::record("Gint::hRGint",this->hRGint->get_memory_size());
		// initialize DMRGint with hRGint when NSPIN != 4
		for (int is = 0; is < this->DMRGint.size(); is++)
		{
			if (this->DMRGint[is] != nullptr)
			{
				delete this->DMRGint[is];
			}
			this->DMRGint[is] = new hamilt::HContainer<double>(*this->hRGint);
		}
		ModuleBase::Memory::record("Gint::DMRGint",this->DMRGint[0]->get_memory_size() * this->DMRGint.size());
	}
	else
	{
		this->hRGintCd->allocate(0);
		ModuleBase::Memory::record("Gint::hRGintCd",this->hRGintCd->get_memory_size());
		for (int is = 0; is < this->DMRGint.size(); is++)
		{
			this->DMRGint[is]->allocate(0);
		}
		ModuleBase::Memory::record("Gint::DMRGint",this->DMRGint[0]->get_memory_size() * this->DMRGint.size());
#ifdef __MPI	
		this->DMRGint_full->allocate(0);
		ModuleBase::Memory::record("Gint::DMRGint_full",this->DMRGint_full->get_memory_size());
#endif
	}

}

void Gint::transfer_DM2DtoGrid(std::vector<hamilt::HContainer<double>*> DM2D)
{
    ModuleBase::TITLE("Gint","transfer_DMR");
    ModuleBase::timer::tick("Gint","transfer_DMR");
	if(GlobalV::NSPIN != 4)
	{
		for (int is = 0; is < this->DMRGint.size(); is++)
		{
#ifdef __MPI
        	hamilt::transferParallels2Serials(*DM2D[is], DMRGint[is]);
#else
			this->DMRGint[is]->set_zero();
			this->DMRGint[is]->add(*DM2D[is]);
#endif	
		}
	}
	else // NSPIN=4 case
	{
#ifdef __MPI
		hamilt::transferParallels2Serials(*DM2D[0], this->DMRGint_full);
#else
		this->DMRGint_full = DM2D[0];
#endif	
		std::vector<double*> tmp_pointer(4, nullptr);
		for(int iap = 0;iap<this->DMRGint_full->size_atom_pairs();++iap)
		{
			auto& ap = this->DMRGint_full->get_atom_pair(iap);
			int iat1 = ap.get_atom_i();
			int iat2 = ap.get_atom_j();
			for(int ir = 0;ir<ap.get_R_size();++ir)
			{
				int* r_index = ap.get_R_index(ir);
				for (int is = 0; is < 4; is++)
				{
					tmp_pointer[is] = this->DMRGint[is]->find_matrix(iat1, iat2, r_index[0], r_index[1], r_index[2])->get_pointer();
				}
				double* data_full = ap.get_pointer(ir);
				for(int irow=0;irow<ap.get_row_size();irow += 2)
				{
					for(int icol=0;icol<ap.get_col_size();icol += 2)
					{
						*(tmp_pointer[0])++ = data_full[icol];
						*(tmp_pointer[1])++ = data_full[icol+1];
					}
					data_full += ap.get_col_size();
					for(int icol=0;icol<ap.get_col_size();icol += 2)
					{
						*(tmp_pointer[2])++ = data_full[icol];
						*(tmp_pointer[3])++ = data_full[icol+1];
					}
					data_full += ap.get_col_size();
				}
			}
		}
	}
    ModuleBase::timer::tick("Gint","transfer_DMR");
}