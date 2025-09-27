#include "symmetry_rho.h"

#include "source_hamilt/module_xc/xc_functional.h"

Symmetry_rho::Symmetry_rho()
{
}

Symmetry_rho::~Symmetry_rho()
{
}

void Symmetry_rho::begin(const int& spin_now,
                         const Charge& chr,
                         const ModulePW::PW_Basis* rho_basis,
                         ModuleSymmetry::Symmetry& symm) const
{
    assert(spin_now < 4); // added by zhengdy-soc

	if (ModuleSymmetry::Symmetry::symm_flag != 1) 
	{
		return;
	}

    ModuleBase::TITLE("Symmetry_rho", "begin");
    ModuleBase::timer::tick("Symmetry_rho","begin");

// both parallel and serial
// if(symm.nrot==symm.nrotk) //pure point-group, do rho_symm in real space
// {
// 	psymm(chr.rho[spin_now], rho_basis, Pgrid, symm);
// 	if(XC_Functional::get_ked_flag()) psymm(chr.kin_r[spin_now],
// rho_basis,Pgrid,symm);
// }
// else	//space group, do rho_symm in reciprocal space

    rho_basis->real2recip(chr.rho[spin_now], chr.rhog[spin_now]);

    psymmg(chr.rhog[spin_now], rho_basis, symm); // need to modify

    rho_basis->recip2real(chr.rhog[spin_now], chr.rho[spin_now]);

    if (XC_Functional::get_ked_flag() || chr.cal_elf)
    {
        // Use std::vector to manage kin_g instead of raw pointer
        std::vector<std::complex<double>> kin_g(chr.ngmc);
        rho_basis->real2recip(chr.kin_r[spin_now], kin_g.data());
        psymmg(kin_g.data(), rho_basis, symm);
        rho_basis->recip2real(kin_g.data(), chr.kin_r[spin_now]);
    }

    ModuleBase::timer::tick("Symmetry_rho","begin");
    return;
}

void Symmetry_rho::begin(const int& spin_now,
                         double** rho,
                         std::complex<double>** rhog,
                         int ngmc,
                         double** kin_r,
                         const ModulePW::PW_Basis* rho_basis,
                         ModuleSymmetry::Symmetry& symm) const
{
    assert(spin_now < 4); // added by zhengdy-soc

    if (ModuleSymmetry::Symmetry::symm_flag != 1)
    {
        return;
    }

    ModuleBase::TITLE("Symmetry_rho", "begin");
    ModuleBase::timer::tick("Symmetry_rho","begin");

    // both parallel and serial
    // if(symm.nrot==symm.nrotk) //pure point-group, do rho_symm in real space
    // {
    // 	psymm(CHR.rho[spin_now], rho_basis, Pgrid, symm);
    // 	if(XC_Functional::get_ked_flag()) psymm(CHR.kin_r[spin_now],
    // rho_basis,Pgrid,symm);
    // }
    // else	//space group, do rho_symm in reciprocal space
    {
        rho_basis->real2recip(rho[spin_now], rhog[spin_now]);
        psymmg(rhog[spin_now], rho_basis, symm);
        rho_basis->recip2real(rhog[spin_now], rho[spin_now]);

        if (XC_Functional::get_ked_flag())
        {
            // Use std::vector to manage kin_g instead of raw pointer
            std::vector<std::complex<double>> kin_g(ngmc);
            rho_basis->real2recip(kin_r[spin_now], kin_g.data());
            psymmg(kin_g.data(), rho_basis, symm);
            rho_basis->recip2real(kin_g.data(), kin_r[spin_now]);
        }
    }

    ModuleBase::timer::tick("Symmetry_rho","begin");
    return;
}

void Symmetry_rho::psymm(double* rho_part,
                         const ModulePW::PW_Basis* rho_basis,
                         Parallel_Grid& Pgrid,
                         ModuleSymmetry::Symmetry& symm) const
{
    ModuleBase::TITLE("Symmetry_rho", "psymm");
    ModuleBase::timer::tick("Symmetry_rho","psymm");

#ifdef __MPI
    // reduce all rho from the first pool.
    std::vector<double> rhotot;
    if (GlobalV::MY_RANK == 0)
    {
        rhotot.resize(rho_basis->nxyz);
        ModuleBase::GlobalFunc::ZEROS(rhotot.data(), rho_basis->nxyz);
    }
    Pgrid.reduce(rhotot.data(), rho_part, false);

    if (GlobalV::MY_RANK == 0)
    {
        symm.rho_symmetry(rhotot.data(), rho_basis->nx, rho_basis->ny, rho_basis->nz);
#else
    symm.rho_symmetry(rho_part, rho_basis->nx, rho_basis->ny, rho_basis->nz);
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

	Pgrid.bcast(rhotot.data(), rho_part, GlobalV::MY_RANK);
#endif

    ModuleBase::timer::tick("Symmetry_rho","psymm");
    return;
}
