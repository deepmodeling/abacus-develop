#include "symm_rho.h"
#include "source_base/parallel_reduce.h"
#include "source_base/parallel_global.h"
#include "source_cell/module_symmetry/symm_rot_spin.h"
#include "source_hamilt/module_xc/xc_functional.h"


void Symmetry_rho::psymmg(std::complex<double>* rhog_part, const ModulePW::PW_Basis *rho_basis, ModuleSymmetry::Symmetry &symm) const
{		
	//(1) get fftixy2is and do Allreduce
	int * fftixy2is = new int [rho_basis->fftnxy];
	rho_basis->getfftixy2is(fftixy2is);		//current proc
#ifdef __MPI
	Parallel_Reduce::reduce_pool(fftixy2is, rho_basis->fftnxy);
	if(rho_basis->poolnproc>1)
		for (int i=0;i<rho_basis->fftnxy;++i)
			fftixy2is[i]+=rho_basis->poolnproc-1;

	// (2) reduce all rho from the first pool.
	std::complex<double>* rhogtot;
	int* ig2isztot = nullptr;
	if(GlobalV::RANK_IN_POOL == 0)
	{
		rhogtot = new std::complex<double>[rho_basis->npwtot];
		ModuleBase::GlobalFunc::ZEROS(rhogtot, rho_basis->npwtot);
		ig2isztot = new int[rho_basis->npwtot];
		ModuleBase::GlobalFunc::ZEROS(rhogtot, rho_basis->npwtot);
	}
	// find max_npw
	int max_npw=0;
	for (int proc = 0; proc < rho_basis->poolnproc; ++proc)
	{
		if(rho_basis->npw_per[proc] > max_npw)
		{
			max_npw=rho_basis->npw_per[proc];
		}
	}
	this->reduce_to_fullrhog(rho_basis, rhogtot, rhog_part, ig2isztot, rho_basis->ig2isz, max_npw);

	// (3) get ixy2ipw and do rhog_symmetry on proc 0 of each pool
	if(GlobalV::RANK_IN_POOL==0)
	{
#endif
		//init ixyz2ipw
		int* ixyz2ipw = new int[rho_basis->fftnxyz];
		for(int i=0;i<rho_basis->fftnxyz;++i) ixyz2ipw[i]=-1;
		// The density must be symmetrized with the same group used to fold the k-points. For
		// nspin=4 magnetic that is the Shubnikov group; Theta leaves the charge invariant, so the
		// antiunitary elements act on rho exactly like unitary ones (their trs_inv is not used here).
		std::vector<ModuleBase::Matrix3> kgmat;
		std::vector<ModuleBase::Vector3<double>> gtr;
		std::vector<double> trs_inv;
		const int nop = symm.density_sym_ops(kgmat, gtr, trs_inv);
#ifdef __MPI
		this->get_ixyz2ipw(rho_basis, ig2isztot, fftixy2is, ixyz2ipw);	
		symm.rhog_symmetry(rhogtot, ixyz2ipw, rho_basis->nx, rho_basis->ny, rho_basis->nz, 
			rho_basis->fftnx, rho_basis->fftny, rho_basis->fftnz,
			rho_basis->gamma_only, kgmat.data(), gtr.data(), nop);
#else
		this->get_ixyz2ipw(rho_basis, rho_basis->ig2isz, fftixy2is, ixyz2ipw);
		symm.rhog_symmetry(rhog_part, ixyz2ipw, rho_basis->nx, rho_basis->ny, rho_basis->nz,
			rho_basis->fftnx, rho_basis->fftny, rho_basis->fftnz,
			rho_basis->gamma_only, kgmat.data(), gtr.data(), nop);
#endif
		delete[] ixyz2ipw;
#ifdef __MPI
	}

	// (4) send the result to other procs in the same pool
	this->rhog_piece_to_all(rho_basis, rhogtot, rhog_part);

	if(GlobalV::RANK_IN_POOL==0)		
	{
		delete[] rhogtot;
		delete[] ig2isztot;
	}
#endif
	delete[] fftixy2is;
	return;
}

void Symmetry_rho::psymmg_soc(std::complex<double>* rhog_x, std::complex<double>* rhog_y,
	std::complex<double>* rhog_z, const ModulePW::PW_Basis* rho_basis, ModuleSymmetry::Symmetry& symm) const
{
	// build the per-operation spin-rotation matrices W(g) from the cartesian rotation
	// gmatc(g) = direct_to_cartesian(gmatrix(g)) = latvec^-1 * gmatrix(g) * latvec.
	auto build_wspin = [&rho_basis, &symm]() {
		const ModuleBase::Matrix3 latvec = rho_basis->latvec;
		const ModuleBase::Matrix3 ilatvec = latvec.Inverse();
		// index [0,nrotk) unitary, [nrotk, nrotk+nrotk_anti) the spatial parts of the
		// antiunitary elements Theta*g -- same layout as density_sym_ops().
		const int na = symm.magnetic_nspin4 ? symm.nrotk_anti : 0;
		// PR 7664 flipped the spinor->Pauli rho^y sign (func_xyz_to_updown), so the grid
		// magnetization rho^y = chr.rho[2] consumed by rhog_symmetry_nspin4 now uses the
		// standard sigma_y=[[0,-i],[i,0]] convention. Relative to that convention the
		// pseudovector rotation applied to (rho^x,rho^y,rho^z) must be conjugated by
		// S=diag(1,-1,1) (i.e. the y-channel handedness is flipped); without this the
		// density symmetrization is out of sync only in the transverse y-channel.
		// spin_so3() itself is left untouched (it is still the physical pseudovector 
		// rotation used for magnetic-group detection on the STRU moments). 
		auto yflip = [](ModuleBase::Matrix3 W) {
			W.e12 = -W.e12; W.e21 = -W.e21; W.e23 = -W.e23; W.e32 = -W.e32; return W; };
		std::vector<ModuleBase::Matrix3> wspin(symm.nrotk + na);
		for (int i = 0; i < symm.nrotk; ++i)
		{
			const ModuleBase::Matrix3 gmatc = ilatvec * symm.gmatrix[i] * latvec;
			wspin[i] = yflip(ModuleSymmetry::SpinRotation::spin_so3(gmatc));
		}
		for (int j = 0; j < na; ++j)
		{
			const ModuleBase::Matrix3 gmatc = ilatvec * symm.gmatrix_anti[j] * latvec;
			wspin[symm.nrotk + j] = yflip(ModuleSymmetry::SpinRotation::spin_so3(gmatc));
		}
		return wspin;
	};

	//(1) get fftixy2is and do Allreduce
	int * fftixy2is = new int [rho_basis->fftnxy];
	rho_basis->getfftixy2is(fftixy2is);		//current proc
#ifdef __MPI
	Parallel_Reduce::reduce_pool(fftixy2is, rho_basis->fftnxy);
	if(rho_basis->poolnproc>1)
		for (int i=0;i<rho_basis->fftnxy;++i)
			fftixy2is[i]+=rho_basis->poolnproc-1;

	// (2) reduce all three spin components from the first pool.
	std::complex<double>* rhogtot_x = nullptr;
	std::complex<double>* rhogtot_y = nullptr;
	std::complex<double>* rhogtot_z = nullptr;
	int* ig2isztot = nullptr;
	if(GlobalV::RANK_IN_POOL == 0)
	{
		rhogtot_x = new std::complex<double>[rho_basis->npwtot];
		rhogtot_y = new std::complex<double>[rho_basis->npwtot];
		rhogtot_z = new std::complex<double>[rho_basis->npwtot];
		ModuleBase::GlobalFunc::ZEROS(rhogtot_x, rho_basis->npwtot);
		ModuleBase::GlobalFunc::ZEROS(rhogtot_y, rho_basis->npwtot);
		ModuleBase::GlobalFunc::ZEROS(rhogtot_z, rho_basis->npwtot);
		ig2isztot = new int[rho_basis->npwtot];
		ModuleBase::GlobalFunc::ZEROS(ig2isztot, rho_basis->npwtot);
	}
	// find max_npw
	int max_npw=0;
	for (int proc = 0; proc < rho_basis->poolnproc; ++proc)
	{
		if(rho_basis->npw_per[proc] > max_npw)
		{
			max_npw=rho_basis->npw_per[proc];
		}
	}
	this->reduce_to_fullrhog(rho_basis, rhogtot_x, rhog_x, ig2isztot, rho_basis->ig2isz, max_npw);
	this->reduce_to_fullrhog(rho_basis, rhogtot_y, rhog_y, ig2isztot, rho_basis->ig2isz, max_npw);
	this->reduce_to_fullrhog(rho_basis, rhogtot_z, rhog_z, ig2isztot, rho_basis->ig2isz, max_npw);

	// (3) get ixy2ipw and do rhog_symmetry_nspin4 on proc 0 of each pool
	if(GlobalV::RANK_IN_POOL==0)
	{
#endif
		//init ixyz2ipw
		int* ixyz2ipw = new int[rho_basis->fftnxyz];
		for(int i=0;i<rho_basis->fftnxyz;++i) ixyz2ipw[i]=-1;
		std::vector<ModuleBase::Matrix3> wspin = build_wspin();
		std::vector<ModuleBase::Matrix3> kgmat;
		std::vector<ModuleBase::Vector3<double>> gtr;
		std::vector<double> trs_inv;
		const int nop = symm.density_sym_ops(kgmat, gtr, trs_inv);
#ifdef __MPI
		this->get_ixyz2ipw(rho_basis, ig2isztot, fftixy2is, ixyz2ipw);
		symm.rhog_symmetry_nspin4(rhogtot_x, rhogtot_y, rhogtot_z, wspin.data(), ixyz2ipw,
			rho_basis->nx, rho_basis->ny, rho_basis->nz,
			rho_basis->fftnx, rho_basis->fftny, rho_basis->fftnz,
			trs_inv.data(), kgmat.data(), gtr.data(), nop);
#else
		this->get_ixyz2ipw(rho_basis, rho_basis->ig2isz, fftixy2is, ixyz2ipw);
		symm.rhog_symmetry_nspin4(rhog_x, rhog_y, rhog_z, wspin.data(), ixyz2ipw,
			rho_basis->nx, rho_basis->ny, rho_basis->nz,
			rho_basis->fftnx, rho_basis->fftny, rho_basis->fftnz,
			trs_inv.data(), kgmat.data(), gtr.data(), nop);
#endif
		delete[] ixyz2ipw;
#ifdef __MPI
	}

	// (4) send the result to other procs in the same pool
	this->rhog_piece_to_all(rho_basis, rhogtot_x, rhog_x);
	this->rhog_piece_to_all(rho_basis, rhogtot_y, rhog_y);
	this->rhog_piece_to_all(rho_basis, rhogtot_z, rhog_z);

	if(GlobalV::RANK_IN_POOL==0)
	{
		delete[] rhogtot_x;
		delete[] rhogtot_y;
		delete[] rhogtot_z;
		delete[] ig2isztot;
	}
#endif
	delete[] fftixy2is;
	return;
}

#ifdef __MPI

void Symmetry_rho::reduce_to_fullrhog(const ModulePW::PW_Basis *rho_basis, 
	std::complex<double>* rhogtot, std::complex<double>* rhogin, 
	int* ig2isztot, const int* ig2iszin, int max_npw) const
{
	ModuleBase::TITLE("Symmetry_rho","reduce_to_fullrhog");

	std::complex<double>* rhog_piece = new std::complex<double>[max_npw];
	int* ig2isz_piece = new int[max_npw];
	
	int npw_start=0;
	for(int proc=0; proc<rho_basis->poolnproc; ++proc)
	{
		ModuleBase::GlobalFunc::ZEROS(rhog_piece, max_npw);
		ModuleBase::GlobalFunc::ZEROS(ig2isz_piece, max_npw);
		
		MPI_Status ierror;

		// case 1: the first part of rho in processor 0 in each pool.
		if(proc == 0 && GlobalV::RANK_IN_POOL ==0)
		{
			for(int ig=0; ig<rho_basis->npw; ++ig)
			{
				rhog_piece[ig] = rhogin[ig];
				ig2isz_piece[ig]=ig2iszin[ig];
			}
		}

		// case 2: > first part rho: send the rho to
		// processor 0 in each pool
		else if(proc == GlobalV::RANK_IN_POOL )
		{
			for(int ig=0; ig<rho_basis->npw; ++ig)
			{
				rhog_piece[ig] = rhogin[ig];
				ig2isz_piece[ig]=ig2iszin[ig];
			}
			MPI_Send(rhog_piece,rho_basis->npw, MPI_DOUBLE_COMPLEX, 0, proc, POOL_WORLD);
			MPI_Send(ig2isz_piece, rho_basis->npw, MPI_INT, 0, proc+rho_basis->poolnproc, POOL_WORLD);
		}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
		else if(GlobalV::RANK_IN_POOL==0)
		{
			MPI_Recv(rhog_piece, rho_basis->npw_per[proc], MPI_DOUBLE_COMPLEX, proc, proc, POOL_WORLD, &ierror);
			MPI_Recv(ig2isz_piece, rho_basis->npw_per[proc], MPI_INT, proc, proc+rho_basis->poolnproc, POOL_WORLD,  &ierror);
		}

		if(GlobalV::RANK_IN_POOL==0)
		{
			for(int ig=0; ig<rho_basis->npw_per[proc]; ++ig)
			{
				rhogtot[npw_start+ig] = rhog_piece[ig];
				ig2isztot[npw_start+ig] = ig2isz_piece[ig];
			}	
			npw_start+=rho_basis->npw_per[proc]; 
		}
	}
	if(GlobalV::RANK_IN_POOL==0)	assert(npw_start==rho_basis->npwtot);
	delete[] rhog_piece;	
	delete[] ig2isz_piece;

	MPI_Barrier(MPI_COMM_WORLD);

	return;
}

void Symmetry_rho::rhog_piece_to_all(const ModulePW::PW_Basis *rho_basis, 
	std::complex<double>* rhogtot, std::complex<double>* rhog_part) const
{	
	ModuleBase::TITLE(" Symmetry_rho","rhog_piece_to_all");

	MPI_Status ierror;

	if(GlobalV::RANK_IN_POOL==0)
	{
		// proc 0: send to other proc in pool
		// itself: directly copy
		for(int ig=0;ig<rho_basis->npw;++ig)
		{
			rhog_part[ig]=rhogtot[ig];
		}
		int npw_start=rho_basis->npw;
		for(int proc=1;proc<rho_basis->poolnproc;++proc)
		{
			MPI_Send(&rhogtot[npw_start], rho_basis->npw_per[proc], MPI_DOUBLE_COMPLEX, proc, proc, POOL_WORLD);
			npw_start+=rho_basis->npw_per[proc];
		}
		assert(npw_start==rho_basis->npwtot);
	}// GlobalV::RANK_IN_POOL == 0
	else
	{
		MPI_Recv(rhog_part, rho_basis->npw_per[GlobalV::RANK_IN_POOL], MPI_DOUBLE_COMPLEX, 0, GlobalV::RANK_IN_POOL, POOL_WORLD, &ierror);
	}
	return;	
}

#endif

// only for MYRANK==0
void Symmetry_rho::get_ixyz2ipw(const ModulePW::PW_Basis *rho_basis, 
	const int* ig2isztot, const int* fftixy2is, int* ixyz2ipw) const
{
	//step 1: get ipsz2ipw
	
	//get ipsz2ipw from ig2isztot
	int* ipsz2ipw = new int [rho_basis->nstot*rho_basis->nz];
	for(int i=0;i<rho_basis->nstot*rho_basis->nz;++i) ipsz2ipw[i]=-1;

	int npw_count=0;
	int nstnz_count=0;
	int ipsz=0;	//global index of a z-grid on stick
	int isz=0;	//local index of a z-grid stick on ip core
	int ipw=0;	// global index of pw (in npwtot)
	for (int ip=0;ip<rho_basis->poolnproc;++ip)
	{
		for (int ig=0;ig<rho_basis->npw_per[ip];++ig)
		{
			ipw=npw_count+ig;
			isz=ig2isztot[ipw];
			ipsz=nstnz_count+isz;
			ipsz2ipw[ipsz]=ipw;
		}
		npw_count+=rho_basis->npw_per[ip];
		nstnz_count+=rho_basis->nst_per[ip]*rho_basis->nz;
	}
	assert(npw_count==rho_basis->npwtot);
	assert(nstnz_count==rho_basis->nstot*rho_basis->nz);

	//step2: ixyz to ipsz

	//save the start-index of (nst*nz) till each core
    int* nstnz_start = new int[rho_basis->poolnproc];
    nstnz_start[0]=0;
    for (int ip=1; ip<rho_basis->poolnproc; ++ip)
        nstnz_start[ip]=nstnz_start[ip-1]+rho_basis->nst_per[ip-1]*rho_basis->nz;

    //tmp variables
    int ixy, ixyz, ip, is, ig=0;
	
	for (int ix=0;ix<rho_basis->fftnx;++ix)
	{
		for (int iy=0;iy<rho_basis->fftny;++iy)
		{
			for(int iz=0;iz<rho_basis->fftnz;++iz)
			{
				ixy = ix*rho_basis->fftny + iy;
				ixyz = ixy*rho_basis->fftnz+iz;
				ip = rho_basis->fftixy2ip[ixy];
				if (ip==-1) continue; //not in any core
				is = fftixy2is[ixy];     //stick-index on ip=proc core
				if (is==-1) continue; //not on any stick
				ipsz = nstnz_start[ip]+is*rho_basis->nz+iz;
				ipw = ipsz2ipw[ipsz];
				ixyz2ipw[ixyz] = ipw;
			}
		}
	}
	assert (ixyz==rho_basis->fftnxyz-1);

	delete[] nstnz_start;
	delete[] ipsz2ipw;
	return;
}
