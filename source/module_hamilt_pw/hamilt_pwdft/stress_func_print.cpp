#include "stress_func.h"
#include "module_base/constants.h"

template<typename FPTYPE, typename Device>
FPTYPE Stress_Func<FPTYPE, Device>::output_acc = 1.0e-8;

//print target stress term 
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::print_stress(const std::string &name, const ModuleBase::matrix& f, const bool screen, bool ry)const
{
	GlobalV::ofs_running << " --------------------------- " << name << " ----------------------------" << std::endl;
	
	
	FPTYPE fac = 1.0;
	
	if(!ry)
	{
	 //	fac = Ry_to_eV / 0.529177;
	}

	std::cout << std::setprecision(5);
	std::cout << std::setiosflags(std::ios::showpos);

	if(screen)
	{
		std::cout << " ------------------- " << name << " --------------------" << std::endl;
	
	}

	for (int i=0;i<3;i++)
	{
		GlobalV::ofs_running << std::setw(15)<< " ";
		if( std::abs(f(i,0)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(i,0) * fac;
		else GlobalV::ofs_running << std::setw(15) << "0";
		if( std::abs(f(i,1)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(i,1) * fac;
		else GlobalV::ofs_running << std::setw(15) << "0";
		if( std::abs(f(i,2)) >output_acc) GlobalV::ofs_running << std::setw(15) << f(i,2) * fac;
		else GlobalV::ofs_running << std::setw(15) << "0";
		GlobalV::ofs_running << std::endl;

		if(screen)
		{
			if( std::abs(f(i,0)) >output_acc) std::cout << std::setw(15) << f(i,0)*fac;
			else std::cout << std::setw(15) << "0";
			if( std::abs(f(i,1)) >output_acc) std::cout << std::setw(15) << f(i,1)*fac;
			else std::cout << std::setw(15) << "0";
			if( std::abs(f(i,2)) >output_acc) std::cout << std::setw(15) << f(i,2)*fac;
			else std::cout << std::setw(15) << "0";
			std::cout << std::endl;
		}	
	}


	std::cout << std::resetiosflags(std::ios::showpos);

    return;
}

//print total stress
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::printstress_total(const ModuleBase::matrix& scs, bool ry)
{
// zhengdy update 2016-10-08
	FPTYPE unit_transform = 1;

	if(!ry)
	{
		unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
	}
//	std::cout.setf(ios::fixed);

	//GlobalV::ofs_running << std::setiosflags(ios::right);
 	GlobalV::ofs_running << std::setprecision(6) << std::setiosflags(std::ios::showpos) << std::setiosflags(std::ios::fixed) << std::endl;
	ModuleBase::GlobalFunc::NEW_PART("TOTAL-STRESS (KBAR)");//Ryd/(a.u.)^3
    std::cout << " ><><><><><><><><><><><><><><><><><><><><><><" << std::endl;
    std::cout << " TOTAL-STRESS (KBAR):" << std::endl;
    std::cout << " ><><><><><><><><><><><><><><><><><><><><><><" << std::endl;

//        if(INPUT.stress_set == 1)
//        int GlobalV::TEST_STRESS = 1;

    
	for (int i=0; i<3; i++)
	{

		//if(GlobalV::TEST_STRESS)
		std::cout << " " << std::setw(15) << scs(i,0)*unit_transform << std::setw(15)
			<< scs(i,1)*unit_transform << std::setw(15) << scs(i,2)*unit_transform << std::endl;

		GlobalV::ofs_running << " " << std::setw(23) << scs(i,0)*unit_transform << std::setw(23)
			<< scs(i,1)*unit_transform << std::setw(23) << scs(i,2)*unit_transform << std::endl;

	}
	FPTYPE pressure = (scs(0,0)+scs(1,1)+scs(2,2))/3.0*unit_transform;
    std::cout << " TOTAL-PRESSURE: " <<pressure<<" KBAR"<< std::endl;
	GlobalV::ofs_running << " TOTAL-PRESSURE: " <<pressure<<" KBAR"<< std::endl;
	GlobalV::ofs_running << std::setiosflags(std::ios::left);
	std::cout << std::resetiosflags(std::ios::showpos);

    return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif