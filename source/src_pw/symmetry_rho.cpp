#include "symmetry_rho.h"
#include "../module_xc/xc_functional.h"
//#include "../src_pw/global.h"

Symmetry_rho::Symmetry_rho()
{

}

Symmetry_rho::~Symmetry_rho()
{

}

void Symmetry_rho::begin(const int &spin_now, const Charge_Broyden &CHR, const PW_Basis &pw, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const
{
	assert(spin_now < 4);//added by zhengdy-soc

	if(!ModuleSymmetry::Symmetry::symm_flag) return;
#ifdef __MPI
	// parallel version
	psymm(CHR.rho[spin_now], pw, Pgrid, symm);
	if(XC_Functional::get_func_type() == 3) psymm(CHR.kin_r[spin_now],pw,Pgrid,symm);
#else
	// series version.
	symm.rho_symmetry(CHR.rho[spin_now], pw.ncx, pw.ncy, pw.ncz);
	if(XC_Functional::get_func_type() == 3) symm.rho_symmetry(CHR.kin_r[spin_now],pw.ncx, pw.ncy, pw.ncz);
#endif
	return;
}

void Symmetry_rho::psymm(double* rho_part, const PW_Basis &pw, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const
{
#ifdef __MPI
	// (1) reduce all rho from the first pool.
	double* rhotot;
	if(GlobalV::MY_RANK == 0)
	{
		rhotot = new double[pw.ncxyz];
		ModuleBase::GlobalFunc::ZEROS(rhotot, pw.ncxyz);
	}
	Pgrid.reduce_to_fullrho(rhotot, rho_part);

	// (2)
	if(GlobalV::MY_RANK==0)
	{
		symm.rho_symmetry(rhotot, pw.ncx, pw.ncy, pw.ncz);
		/*
		int count = 0;
		GlobalV::ofs_running << scientific;
		for(int iz=0; iz<pw.ncz; iz++)
		{
			GlobalV::ofs_running << "\n iz=" << iz;
			for(int iy=0; iy<pw.ncy; iy++)
			{
				for(int ix=0; ix<pw.ncx; ix++)
				{
					if(count%5==0) GlobalV::ofs_running << "\n";
					++count;
					GlobalV::ofs_running << " " << rhotot[ix*pw.ncy*pw.ncz+iy*pw.ncz+iz];
				}
			}
		}
		*/
	}

	// (3)
	const int ncxy = pw.ncx * pw.ncy;
	double* zpiece = new double[ncxy];
	for(int iz=0; iz<pw.ncz; iz++)
	{
		//GlobalV::ofs_running << "\n iz=" << iz;
		ModuleBase::GlobalFunc::ZEROS(zpiece, ncxy);
		if(GlobalV::MY_RANK==0)
		{
			for(int ix=0; ix<pw.ncx; ix++)
			{
				for(int iy=0; iy<pw.ncy; iy++)
				{
					const int ir = ix * pw.ncy + iy;
					zpiece[ir] = rhotot[ix * pw.ncy * pw.ncz + iy * pw.ncz + iz];
					//rho[ir*nczp+znow] = zpiece[ir];
				}
			}
		}
		Pgrid.zpiece_to_all(zpiece,iz, rho_part);
	}

	if(GlobalV::MY_RANK==0)		delete[] rhotot;
	delete[] zpiece;
#endif
	return;
}
