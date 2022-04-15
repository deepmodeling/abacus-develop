#include "istate_envelope.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "src_io/wf_io.h"
#include "src_io/write_wfc_realspace.h"

IState_Envelope::IState_Envelope()
{}

IState_Envelope::~IState_Envelope()
{}


void IState_Envelope::begin(Local_Orbital_wfc &lowf, Gint_Gamma &gg, int& out_wfc_pw, int& out_wfc_r)
{
	ModuleBase::TITLE("IState_Envelope","begin");

	std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

	if(!GlobalV::GAMMA_ONLY_LOCAL)
	{
		ModuleBase::WARNING_QUIT("IState_Envelope::begin","Only available for GAMMA_ONLY_LOCAL now.");
	}

	// (1) 
	// (1.1) allocate the space for GlobalC::LOWF.WFC_GAMMA

	// (1.2) read in LOWF_GAMMA.dat

	// mohan update 2011-03-21
	// if ucell is odd, it's correct,
	// if ucell is even, it's also correct.
	// +1.0e-8 in case like (2.999999999+1)/2
	int fermi_band = static_cast<int>( (GlobalC::CHR.nelec+1)/2 + 1.0e-8 ) ;
	int bands_below = GlobalV::NBANDS_ISTATE;
	int bands_above = GlobalV::NBANDS_ISTATE;

	std::cout << " number of electrons = " << GlobalC::CHR.nelec << std::endl;
	std::cout << " number of occupied bands = " << fermi_band << std::endl;
	std::cout << " plot band decomposed charge density below fermi surface with " 
	<< bands_below << " bands." << std::endl;
	
	std::cout << " plot band decomposed charge density above fermi surface with " 
	<< bands_above << " bands." << std::endl;
	
	// (2) cicle:
	
	// (2.1) calculate the selected density matrix
	// from wave functions.

	// (2.2) carry out the grid integration to
	// get the charge density.

	// (2.3) output the charge density in .cub format.
	this->bands_picked = new bool[GlobalV::NBANDS];
	ModuleBase::GlobalFunc::ZEROS(bands_picked, GlobalV::NBANDS);
	for(int ib=0; ib<GlobalV::NBANDS; ib++)
	{
		if( ib >= fermi_band - bands_below ) 
		{
			if( ib < fermi_band + bands_above)
			{
				bands_picked[ib] = true;
			}
		}
	}

	//allocate grid wavefunction for gamma_only
	std::vector<double**> wfc_gamma_grid(GlobalV::NSPIN);
	for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        wfc_gamma_grid[is] = new double* [GlobalV::NBANDS];
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            wfc_gamma_grid[is][ib] = new double[GlobalC::GridT.lgd];
    }

    //for pw-wfc in G space
    ModuleBase::ComplexMatrix* pw_wfc_g;
    
    if (out_wfc_pw || out_wfc_r)
    {
        pw_wfc_g = new ModuleBase::ComplexMatrix[GlobalC::kv.nks];
        for (int ik = 0;ik < GlobalC::kv.nks;++ik)
            pw_wfc_g[ik].create(GlobalV::NBANDS, GlobalC::kv.ngk[ik], true);
    }


    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
	{
		if(bands_picked[ib])
		{
			for(int is=0; is<GlobalV::NSPIN; ++is)
			{
				std::cout << " Perform envelope function for band " << ib+1 << std::endl;
				ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is],GlobalC::pw.nrxx);	


				//---------------------------------------------------------
				// GlobalC::LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
				// and 2d-to-grid conversion is unified into `wfc_2d_to_grid`.
                //---------------------------------------------------------
#ifdef __MPI
                lowf.wfc_2d_to_grid(0, lowf.wfc_gamma[is].c, wfc_gamma_grid[is]);
#else
                for (int i = 0;i < GlobalV::NBANDS;++i)
                {
                    for (int j = 0;j < GlobalV::NLOCAL;++j)
                        wfc_gamma_grid[is][i][j] = lowf.wfc_gamma[is](i, j);
                }
#endif
				gg.cal_env( wfc_gamma_grid[is][ib], GlobalC::CHR.rho[is] );


				GlobalC::CHR.save_rho_before_sum_band(); //xiaohui add 2014-12-09
				std::stringstream ss;
				ss << GlobalV::global_out_dir << "BAND" << ib + 1 << "_ENV" << is+1 << "_CHG";
				// 0 means definitely output charge density.
				bool for_plot = true;
                GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ss.str(), 3, for_plot);
                
                if (out_wfc_pw || out_wfc_r) //only for gamma_only now
                    this->set_pw_wfc(GlobalC::pw, 0, ib, GlobalV::NSPIN, GlobalC::kv.ngk[0],
                        GlobalC::CHR.rho_save, pw_wfc_g[0]);
            }
		}
	}

    if (out_wfc_pw || out_wfc_r)
    {
        if (out_wfc_pw)
        {
            std::stringstream ssw;
            ssw << GlobalV::global_out_dir << "WAVEFUNC";
            std::cout << " write G-space wavefunction into \"" <<
                GlobalV::global_out_dir << "/" << ssw.str() << "\" files." << std::endl;
            WF_io::write_wfc2(ssw.str(), pw_wfc_g, GlobalC::pw.gcar);
        }
        if(out_wfc_r)
		    Write_Wfc_Realspace::write_wfc_realspace_1(pw_wfc_g, "wfc_realspace", false);
        delete[] pw_wfc_g;
    }  
    
    delete[] bands_picked;
    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        for (int ib = 0;ib < GlobalV::NBANDS; ++ib)
            delete[] wfc_gamma_grid[is][ib];
        delete[] wfc_gamma_grid[is];
    }
    return;
}

//for each band
void IState_Envelope::set_pw_wfc(PW_Basis& pwb,
    const int& ik, const int& ib, const int& nspin,const int& ngk, 
    const double* const* const rho,
    ModuleBase::ComplexMatrix &wfc_g)
{
    if (ib == 0)//once is enough
        ModuleBase::TITLE("IState_Envelope", "set_pw_wfc");
    
    std::vector<std::complex<double>> Porter(pwb.nrxx);
     // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir=0; ir<pwb.nrxx; ir++) 
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);

    //call FFT
    pwb.FFT_wfc.FFT3D(Porter.data(), -1);

    // set pw_wfc_g
    // ig2fftw: the index map from i_ngk(local) to i_ngmw(local)
    for (int ig = 0;ig < ngk;++ig)     // is it right (check local index)???
        wfc_g(ib,ig) = Porter[pwb.ig2fftw[GlobalC::wf.igk(ik, ig)]];
}
