#include "H_Hartree_pw.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"

double H_Hartree_pw::hartree_energy=0.0;

//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
ModuleBase::matrix H_Hartree_pw::v_hartree(
	const UnitCell &cell, 
	PW_Basis &pwb, 
	const int &nspin,
	const double*const*const rho)
{
    ModuleBase::TITLE("H_Hartree_pw","v_hartree");
    ModuleBase::timer::tick("H_Hartree_pw","v_hartree");

    //  Hartree potential VH(r) from n(r)
    std::vector<std::complex<double>> Porter(pwb.nrxx);
    const int nspin0 = (nspin==2) ? 2 : 1;
    for(int is=0; is<nspin0; is++)
        for (int ir=0; ir<pwb.nrxx; ir++) 
            Porter[ir] += std::complex<double>( rho[is][ir], 0.0 );
    //=============================
    //  bring rho (aux) to G space
    //=============================
    pwb.FFT_chg.FFT3D(Porter.data(), -1);

    //double charge;
    //if (pwb.gstart == 1)
    //    charge = cell.omega * Porter[pwb.ig2fftc[0]].real();
    //OUT(GlobalV::ofs_running, "v_h charge", charge);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================

	double ehart = 0.0;

    std::vector<std::complex<double>> vh_g(pwb.ngmc);
    for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * pwb.gg [ig]);

            ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
            vh_g[ig] = fac * Porter[j];
        }
    }

    Parallel_Reduce::reduce_double_pool( ehart );
    ehart *= 0.5 * cell.omega;
    //std::cout << " ehart=" << ehart << std::endl;
    H_Hartree_pw::hartree_energy = ehart;

    std::fill( Porter.begin(), Porter.end(), std::complex<double>(0.0,0.0) );
    for (int ig = 0;ig < pwb.ngmc;ig++)
        Porter[pwb.ig2fftc[ig]] = vh_g[ig];
    //==========================================
    //transform hartree potential to real space
    //==========================================
    pwb.FFT_chg.FFT3D(Porter.data(), 1);

    //==========================================
    //Add hartree potential to the xc potential
    //==========================================
    ModuleBase::matrix v(nspin, pwb.nrxx);
    if(nspin==4)
    {
        for (int ir = 0;ir < pwb.nrxx;ir++)
            v(0, ir) = Porter[ir].real();
    }
    else
    {
        for (int is = 0;is < nspin;is++)
            for (int ir = 0;ir < pwb.nrxx;ir++)
                v(is, ir) = Porter[ir].real();
    }

//-----------------------------------------------------------
// we need to add this out_pot funciton back 
// in near future, 2021-02-25
//-----------------------------------------------------------
	//-------------------------------------------
	// output the Hartree potential into a file.
	//-------------------------------------------
/*
	if(out_pot==-2)
	{
		std::cout << " output VH" << std::endl;
		int is = 0;
		int iter = 0;
		int precision = 3;
		std::string fn = "VH.dat";
		std::stringstream ss;
		ss << GlobalV::global_out_dir << fn;
		ModuleBase::matrix v;
		v.create(1,pwb.nrxx);
		for(int ir=0; ir<pwb.nrxx; ++ir)
		{
			v(0,ir) = Porter[ir].real();
		}
		this->write_potential( is, iter, ss.str(), v, precision, 1 );
	}
*/

    ModuleBase::timer::tick("H_Hartree_pw","v_hartree");
    return v;
} // end subroutine v_h
