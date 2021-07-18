#include "symmetry_rho.h"
//#include "../src_pw/global.h"

Symmetry_rho::Symmetry_rho()
{

}

Symmetry_rho::~Symmetry_rho()
{

}

void Symmetry_rho::begin(const int &spin_now, const Charge_Broyden &CHR, const PW_Basis &pw, const Parallel_Grid &Pgrid, const Symmetry &symm) const
{
	assert(spin_now < 4);//added by zhengdy-soc

	if(!Symmetry::symm_flag) return;
#ifdef __MPI
	// parallel version
	psymm(CHR.rho[spin_now], pw, Pgrid, symm);
#else
	// series version.
	symm.rho_symmetry(CHR.rho[spin_now], pw.ncx, pw.ncy, pw.ncz);
#endif
	return;
}

void Symmetry_rho::psymm(double* rho_part, const PW_Basis &pw, const Parallel_Grid &Pgrid, const Symmetry &symm) const
{
#ifdef __MPI
	// (1) reduce all rho from the first pool.
	double* rhotot = new double[pw.ncxyz];
	ZEROS(rhotot, pw.ncxyz);
	Pgrid.reduce_to_fullrho(rhotot, rho_part);
	
	// (2)
	if(RANK_IN_POOL==0)
	{
		symm.rho_symmetry(rhotot, pw.ncx, pw.ncy, pw.ncz);
		/*
		int count = 0;
		ofs_running << scientific;
		for(int iz=0; iz<pw.ncz; iz++)
		{
			ofs_running << "\n iz=" << iz;
			for(int iy=0; iy<pw.ncy; iy++)
			{
				for(int ix=0; ix<pw.ncx; ix++)
				{
					if(count%5==0) ofs_running << "\n";
					++count;
					ofs_running << " " << rhotot[ix*pw.ncy*pw.ncz+iy*pw.ncz+iz];
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
		//ofs_running << "\n iz=" << iz;
		ZEROS(zpiece, ncxy);
		if(MY_RANK==0)
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
						
	delete[] rhotot;
	delete[] zpiece;
#endif
	return;
}
