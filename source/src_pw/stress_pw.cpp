#include "./stress_pw.h"
#include "vdwd2.h"
#include "vdwd3.h"
#include "../module_base/timer.h"
#include "global.h"

void Stress_PW::cal_stress(ModuleBase::matrix& sigmatot)
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
	stress_kin(sigmakin);
	
	//hartree contribution
	stress_har(sigmahar, 1);

    //ewald contribution
    stress_ewa(sigmaewa, 1);

    //xc contribution: add gradient corrections(non diagonal)
    for(int i=0;i<3;i++)
	{
       sigmaxc(i,i) = - (GlobalC::en.etxc - GlobalC::en.vtxc) / GlobalC::ucell.omega;
    }
    stress_gga(sigmaxc);
    if(XC_Functional::get_func_type() == 3) stress_mgga(sigmaxc);

    //local contribution
    stress_loc(sigmaloc, 1);
    
    //nlcc
    stress_cc(sigmaxcc, 1);
   
    //nonlocal
	stress_nl(sigmanl);

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
    
	if(ModuleSymmetry::Symmetry::symm_flag)                          
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

void Stress_PW::stress_vdw(ModuleBase::matrix& sigma)
{
	ModuleBase::matrix force;
	if(GlobalC::vdwd2_para.flag_vdwd2) //Peize Lin add 2014-04-04, update 2021-03-09
	{
		Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
		vdwd2.cal_stress();
		sigma = vdwd2.get_stress().to_matrix();
	}
	if(GlobalC::vdwd3_para.flag_vdwd3) //jiyy add 2019-05-18, update 2021-05-02
	{
		Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
		vdwd3.cal_stress();
		sigma = vdwd3.get_stress().to_matrix();
	}              
	return;
}
