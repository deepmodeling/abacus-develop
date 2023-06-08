#include "symmetry_rho.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

Symmetry_rho::Symmetry_rho()
{

}

Symmetry_rho::~Symmetry_rho()
{

}

void Symmetry_rho::begin(const int& spin_now,
                         const Charge& CHR,
                         const ModulePW::PW_Basis* rho_basis,
                         Parallel_Grid& Pgrid,
                         ModuleSymmetry::Symmetry& symm) const
{
	assert(spin_now < 4);//added by zhengdy-soc
    if (ModuleSymmetry::Symmetry::symm_flag != 1) return;
    bool &sreal = symm.rhosym_real;
    // both parallel and serial
    if (sreal) //pure point-group, do rho_symm in real space
	{
		sreal = psymm(CHR.rho[spin_now], rho_basis, Pgrid, symm);
        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5 && sreal)
            sreal = psymm(CHR.kin_r[spin_now], rho_basis, Pgrid, symm);
        if (sreal) return;
        else GlobalV::ofs_running << "WARNING: too large error in real-space rho_symmetry, do it in reciprocal space."<<std::endl;
    }
    //space group, or unacceptable error in real-space rho_symm: 
    //do rho_symm in reciprocal space
    rho_basis->real2recip(CHR.rho[spin_now], CHR.rhog[spin_now]);
    psymmg(CHR.rhog[spin_now], rho_basis, Pgrid, symm);	//need to modify
    rho_basis->recip2real(CHR.rhog[spin_now], CHR.rho[spin_now]);

    if(XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) 
    {
        std::complex<double>* kin_g = new std::complex<double>[CHR.ngmc];
        rho_basis->real2recip(CHR.kin_r[spin_now], kin_g);
        psymmg(kin_g,  rho_basis,Pgrid,symm);
        rho_basis->recip2real(kin_g, CHR.kin_r[spin_now]);
        delete[] kin_g;
    }	
	return;
}

bool Symmetry_rho::psymm(double* rho_part, const ModulePW::PW_Basis *rho_basis, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const
{
    bool sreal = true;
#ifdef __MPI
	// (1) reduce all rho from the first pool.
	double* rhotot;
	if(GlobalV::MY_RANK == 0)
	{
		rhotot = new double[rho_basis->nxyz];
		ModuleBase::GlobalFunc::ZEROS(rhotot, rho_basis->nxyz);
	}
	Pgrid.reduce_to_fullrho(rhotot, rho_part);

	// (2)
	if(GlobalV::MY_RANK==0)
	{
		sreal = symm.rho_symmetry(rhotot, rho_basis->nx, rho_basis->ny, rho_basis->nz);
#else
		sreal = symm.rho_symmetry(rho_part, rho_basis->nx, rho_basis->ny, rho_basis->nz);
#endif
		/*
		int count = 0;
		GlobalV::ofs_running << scientific;
		for(int iz=0; iz<rho_basis->nz; iz++)
		{
			GlobalV::ofs_running << "\n iz=" << iz;
			for(int iy=0; iy<rho_basis->ny; iy++)
			{
				for(int ix=0; ix<rho_basis->nx; ix++)
				{
					if(count%5==0) GlobalV::ofs_running << "\n";
					++count;
					GlobalV::ofs_running << " " << rhotot[ix*rho_basis->ny*rho_basis->nz+iy*rho_basis->nz+iz];
				}
			}
		}
		*/
	#ifdef __MPI
	}

	// (3)
	const int ncxy = rho_basis->nx * rho_basis->ny;
	double* zpiece = new double[ncxy];
    if(sreal)
	for(int iz=0; iz<rho_basis->nz; iz++)
	{
		//GlobalV::ofs_running << "\n iz=" << iz;
		ModuleBase::GlobalFunc::ZEROS(zpiece, ncxy);
		if(GlobalV::MY_RANK==0)
		{
			for(int ix=0; ix<rho_basis->nx; ix++)
			{
				for(int iy=0; iy<rho_basis->ny; iy++)
				{
					const int ir = ix * rho_basis->ny + iy;
					zpiece[ir] = rhotot[ix * rho_basis->ny * rho_basis->nz + iy * rho_basis->nz + iz];
					//rho[ir*nczp+znow] = zpiece[ir];
				}
			}
		}
		Pgrid.zpiece_to_all(zpiece,iz, rho_part);
	}

	if(GlobalV::MY_RANK==0)		delete[] rhotot;
	delete[] zpiece;
#endif
	return sreal;
}
