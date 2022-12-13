#include "./stress_pw.h"
#include "../module_base/timer.h"
#include "global.h"
#include "module_vdw/vdw.h"

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::cal_stress(ModuleBase::matrix& sigmatot, const psi::Psi<complex<FPTYPE>>* psi_in, const psi::Psi<complex<FPTYPE>, Device>* d_psi_in)
{
	ModuleBase::TITLE("Stress_PW","cal_stress");
	ModuleBase::timer::tick("Stress_PW","cal_stress");    

	// total stress
	sigmatot.create(3,3);
	ModuleBase::matrix sigmaxc;
	// exchange-correlation stress
	sigmaxc.create(3,3);
	// hartree stress
	ModuleBase::matrix sigmahar;
	sigmahar.create(3,3);
	// electron kinetic stress
	ModuleBase::matrix sigmakin;
	sigmakin.create(3,3);
	// local pseudopotential stress
	ModuleBase::matrix sigmaloc;
	sigmaloc.create(3,3);
	// non-local pseudopotential stress
	ModuleBase::matrix sigmanl;
	sigmanl.create(3,3);
	// Ewald stress
	ModuleBase::matrix sigmaewa;
	sigmaewa.create(3,3);
	// non-linear core correction stress
	ModuleBase::matrix sigmaxcc;
	sigmaxcc.create(3,3);
	// vdw stress
	ModuleBase::matrix sigmavdw;
	sigmavdw.create(3,3);

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			sigmatot(i,j) = 0.0;
			sigmaxc(i,j) = 0.0;
			sigmahar(i,j) = 0.0;
			sigmakin(i,j) = 0.0;
			sigmaloc(i,j) = 0.0;
			sigmanl(i,j) = 0.0;
			sigmaewa(i,j) = 0.0;
			sigmaxcc(i,j) = 0.0;
			sigmavdw(i,j) = 0.0;
		}
	}

	//kinetic contribution
	this->stress_kin(sigmakin, this->pelec->wg, psi_in);
	
	//hartree contribution
    this->stress_har(sigmahar, GlobalC::rhopw, 1, pelec->charge);

    //ewald contribution
    this->stress_ewa(sigmaewa, GlobalC::rhopw, 1);

    //xc contribution: add gradient corrections(non diagonal)
    for(int i=0;i<3;i++)
	{
       sigmaxc(i,i) = - (GlobalC::en.etxc - GlobalC::en.vtxc) / GlobalC::ucell.omega;
    }
    this->stress_gga(sigmaxc, pelec->charge);
    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5) {
        this->stress_mgga(sigmaxc, this->pelec->wg, this->pelec->pot->get_effective_vofk(), pelec->charge, psi_in);
    }

    //local contribution
    this->stress_loc(sigmaloc, GlobalC::rhopw, 1, pelec->charge);
    
    //nlcc
    this->stress_cc(sigmaxcc, GlobalC::rhopw, 1, pelec->charge);
   
    //nonlocal
	this->stress_nl(sigmanl, this->pelec->wg, d_psi_in);

	//vdw term
	stress_vdw(sigmavdw);

    for(int ipol=0;ipol<3;ipol++)
	{
        for(int jpol=0;jpol<3;jpol++)
		{
			sigmatot(ipol,jpol) = sigmakin(ipol,jpol) 
								+ sigmahar(ipol,jpol) 
								+ sigmanl(ipol,jpol) 
								+ sigmaxc(ipol,jpol) 
								+ sigmaxcc(ipol,jpol) 
								+ sigmaewa(ipol,jpol)
								+ sigmaloc(ipol,jpol)
								+ sigmavdw(ipol,jpol);
        }
    }
    
	if(ModuleSymmetry::Symmetry::symm_flag == 1)                          
	{
		GlobalC::symm.stress_symmetry(sigmatot, GlobalC::ucell);
	}

	bool ry = false;
	this->printstress_total(sigmatot, ry);

	if(GlobalV::TEST_STRESS) 
	{               
		GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
		GlobalV::ofs_running << std::setiosflags(ios::showpos);
		GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(10) << std::endl;
		this->print_stress("KINETIC    STRESS",sigmakin,GlobalV::TEST_STRESS,ry);
		this->print_stress("LOCAL    STRESS",sigmaloc,GlobalV::TEST_STRESS,ry);
		this->print_stress("HARTREE    STRESS",sigmahar,GlobalV::TEST_STRESS,ry);
		this->print_stress("NON-LOCAL    STRESS",sigmanl,GlobalV::TEST_STRESS,ry);
		this->print_stress("XC    STRESS",sigmaxc,GlobalV::TEST_STRESS,ry);
		this->print_stress("EWALD    STRESS",sigmaewa,GlobalV::TEST_STRESS,ry);
		this->print_stress("NLCC    STRESS",sigmaxcc,GlobalV::TEST_STRESS,ry);
		this->print_stress("TOTAL    STRESS",sigmatot,GlobalV::TEST_STRESS,ry);
	}
	ModuleBase::timer::tick("Stress_PW","cal_stress");
	return;
    
}

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_vdw(ModuleBase::matrix& sigma)
{
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, INPUT);
    if (vdw_solver != nullptr)
    {
    sigma = vdw_solver->get_stress().to_matrix();
    }
	return;
}

template class Stress_PW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_PW<double, psi::DEVICE_GPU>;
#endif
