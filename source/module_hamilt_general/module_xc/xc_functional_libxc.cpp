#ifdef USE_LIBXC

#include "xc_functional_libxc.h"
#include "module_parameter/parameter.h"
#include "module_base/tool_quit.h"

#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"		// just for GlobalC::exx_info
#endif

#include <xc.h>

#include <vector>

std::pair<int,std::vector<int>> XC_Functional_Libxc::set_xc_type_libxc(std::string xc_func_in)
{
    // determine the type (lda/gga/mgga)
	int func_type; //0:none, 1:lda, 2:gga, 3:mgga, 4:hybrid lda/gga, 5:hybrid mgga
    func_type = 1;
    if(xc_func_in.find("GGA") != std::string::npos) { func_type = 2; }
    if(xc_func_in.find("MGGA") != std::string::npos) { func_type = 3; }
    if(xc_func_in.find("HYB") != std::string::npos) { func_type =4; }
    if(xc_func_in.find("HYB") != std::string::npos && xc_func_in.find("MGGA") != std::string::npos) { func_type =5; }

    // determine the id
	std::vector<int> func_id; // libxc id of functional
    int pos = 0;
    std::string delimiter = "+";
    std::string token;
    while ((pos = xc_func_in.find(delimiter)) != std::string::npos)
    {
        token = xc_func_in.substr(0, pos);
        int id = xc_functional_get_number(token.c_str());
        std::cout << "func,id" << token << " " << id << std::endl;
        if (id == -1) { ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc","functional name not recognized!"); }
        func_id.push_back(id);
        xc_func_in.erase(0, pos + delimiter.length());
    }
    int id = xc_functional_get_number(xc_func_in.c_str());
    std::cout << "func,id" << xc_func_in << " " << id << std::endl;
    if (id == -1) { ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc","functional name not recognized!"); }
    func_id.push_back(id);

    return std::make_pair(func_type, func_id);
}

std::vector<xc_func_type> XC_Functional_Libxc::init_func(const std::vector<int> &func_id, const int xc_polarized)
{
	// 'funcs' is the return value
	std::vector<xc_func_type> funcs;

	//-------------------------------------------
	// define a function named 'add_func', which
	// initialize a functional according to its ID
	//-------------------------------------------
	auto add_func = [&]( const int func_id )
	{
		funcs.push_back({});
		// 'xc_func_init' is defined in Libxc
		xc_func_init( &funcs.back(), func_id, xc_polarized );
	};

	for(int id : func_id)
	{
        if(id == XC_LDA_XC_KSDT || id == XC_LDA_XC_CORRKSDT || id == XC_LDA_XC_GDSMFB) //finite temperature XC functionals
        {
            add_func(id);
            double parameter_finitet[1] = {PARAM.inp.xc_temperature * 0.5}; // converts to Hartree for libxc
            xc_func_set_ext_params(&funcs.back(), parameter_finitet);
        }
#ifdef __EXX
		else if( id == XC_HYB_GGA_XC_PBEH ) // PBE0
		{
			add_func( XC_HYB_GGA_XC_PBEH );
			double parameter_hse[3] = { GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_info.info_global.hse_omega,
				GlobalC::exx_info.info_global.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);
		}
		else if( id == XC_HYB_GGA_XC_HSE06 ) // HSE06 hybrid functional
		{
			add_func( XC_HYB_GGA_XC_HSE06 );
			double parameter_hse[3] = { GlobalC::exx_info.info_global.hybrid_alpha,
				GlobalC::exx_info.info_global.hse_omega,
				GlobalC::exx_info.info_global.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);
		}
        // added by jghan, 2024-07-06
		else if( id == XC_GGA_X_ITYH ) // short-range of B88_X
		{
			add_func( XC_GGA_X_ITYH );
			double parameter_omega[1] = {GlobalC::exx_info.info_global.hse_omega};
			xc_func_set_ext_params(&funcs.back(), parameter_omega);	
		}
		else if( id == XC_GGA_C_LYPR ) // short-range of LYP_C
		{
			add_func( XC_GGA_C_LYPR );
            // the first six parameters come from libxc, and may need to be modified in some cases
			double parameter_lypr[7] = {0.04918, 0.132, 0.2533, 0.349, 0.35/2.29, 2.0/2.29, GlobalC::exx_info.info_global.hse_omega};
			xc_func_set_ext_params(&funcs.back(), parameter_lypr);	
		}
#endif
		else
		{
			add_func( id );
		}
	}
	return funcs;
}

void XC_Functional_Libxc::finish_func(std::vector<xc_func_type> &funcs)
{
    for(xc_func_type func : funcs)
	{
        xc_func_end(&func);
    }
}

#endif