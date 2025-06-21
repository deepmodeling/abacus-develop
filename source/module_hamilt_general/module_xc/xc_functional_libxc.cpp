#ifdef USE_LIBXC

#include "xc_functional_libxc.h"
#include "module_parameter/parameter.h"
#include "source_base/tool_quit.h"
#include "source_base/formatter.h"

#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/global.h"		// just for GlobalC::exx_info
#endif

#include <xc.h>
#include <vector>
#include <regex>
#include <map>
#include <algorithm>
#include <cassert>

bool not_supported_xc_with_laplacian(const std::string& xc_func_in)
{
	// see Pyscf: https://github.com/pyscf/pyscf/blob/master/pyscf/dft/libxc.py#L1062
	// ABACUS issue: https://github.com/deepmodeling/abacus-develop/issues/5372
	const std::vector<std::string> not_supported = {
		"MGGA_XC_CC06", "MGGA_C_CS", "MGGA_X_BR89", "MGGA_X_MK00"};
	for (const std::string& s : not_supported)
	{
		if (xc_func_in.find(s) != std::string::npos)
		{
			return true;
		}
	}
	return false;
}

bool not_supported_xc_with_nonlocal_vdw(const std::string& xc_func_in)
{
	const std::string xc_func = FmtCore::upper(xc_func_in);
	if(xc_func.find("VDW") != std::string::npos) { return true; }
	/* known excluded: GGA_X_OPTB86B_VDW, GGA_X_OPTB88_VDW, GGA_X_OPTPBE_VDW, GGA_X_PBEK1_VDW */

	if(xc_func.find("VV10") != std::string::npos) { return true; }
	/* known excluded: GGA_XC_VV10, HYB_GGA_XC_LC_VV10, MGGA_C_REVSCAN_VV10, MGGA_C_SCAN_VV10, 
	            	   MGGA_C_SCANL_VV10, MGGA_XC_VCML_RVV10 */
					   
	const std::vector<std::string> not_supported = {"C09X", "VCML", "HYB_MGGA_XC_WB97M_V", "MGGA_XC_B97M_V"};
	for(const std::string& str : not_supported)
	{
		if(xc_func.find(str) != std::string::npos) { return true; }
	}
	/* known excluded: GGA_X_C09X, MGGA_X_VCML, HYB_MGGA_XC_WB97M_V, MGGA_XC_B97M_V */

	/* There is also a functional not quite sure: HYB_GGA_XC_WB97X_V */
	if(xc_func.find("HYB_GGA_XC_WB97X_V") != std::string::npos)
	{
		std::cout << " WARNING: range-seperated XC omega-B97 family with nonlocal correction term is used.\n" 
		          << "          if you are not planning to use these functionals like wB97X-D3BJ that:\n"
				  << "          XC_GGA_XC_WB97X_V with specified D3BJ DFT-D3 parameters, this is not what\n"
				  << "          you want." << std::endl;
	}
	return false;
}

int xc_func_type_classifier(const std::string& xc_func, 
							const std::map<std::string, int>& mymap = {
								{"LDA", 1},
								{"GGA", 2},
								{"MGGA", 3},
								{"HYB_LDA", 4},
								{"HYB_GGA", 4},
								{"HYB_MGGA", 5}
							})
{
	// the libxc standard functional pattern is like:
	// "(XC_)?(LDA|GGA|MGGA|HYB_GGA|HYB_MGGA|HYB_LDA)_(X|C|XC|K)(_(.*))?"
	std::regex pattern(R"((XC_)?(LDA|GGA|MGGA|HYB_GGA|HYB_MGGA|HYB_LDA)_(X|C|XC|K)(_(.*))?)");
	std::smatch match;
	if (std::regex_match(xc_func, match, pattern)) {
		std::string type = match[2].str();
		auto it = mymap.find(type);
		if (it != mymap.end()) {
			return it->second;
		} else {
			ModuleBase::WARNING_QUIT("XC_Functional_Libxc::xc_func_type_classifier",
				"Unrecognized functional type: " + type);
		}
	} else {
		ModuleBase::WARNING_QUIT("XC_Functional_Libxc::xc_func_type_classifier",
			"Unrecognized functional format: " + xc_func);
	}
}

std::pair<int,std::vector<int>> 
XC_Functional_Libxc::set_xc_type_libxc(const std::string& xc_func_in)
{
    // check if the functional involves Laplacian of rho
	if (not_supported_xc_with_laplacian(xc_func_in))
	{
		ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc",
			"XC Functional involving Laplacian of rho is not implemented.");
	}

	// check if the functional involves non-local dispersion
	if(not_supported_xc_with_nonlocal_vdw(xc_func_in))
	{ 
		ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc",
			"functionals with non-local dispersion are not supported."); 
	}

	// check the consistency of the functional type (LDA, GGA, MGGA, HYB_LDA, HYB_GGA, HYB_MGGA)
	const std::vector<std::string> xcfunc_words_ = FmtCore::split(xc_func_in, "+");
	std::vector<int> xcfunc_type_(xcfunc_words_.size(), 0); // 0: None
	std::transform(xcfunc_words_.begin(), xcfunc_words_.end(), xcfunc_type_.begin(),
		[](const std::string& func) { return xc_func_type_classifier(func); });
	if (std::adjacent_find(xcfunc_type_.begin(), xcfunc_type_.end(),
		[](int a, int b) { return a != b; }) != xcfunc_type_.end())
	{
		ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc",
			"All exchange-correlation functionals must be of the same type"
			"(LDA, GGA, MGGA, HYB_LDA, HYB_GGA, HYB_MGGA).");
	}

	// check if there is None (no, we dont check it)
    int func_type = xcfunc_type_.front(); // all functionals are of the same type
	// if (func_type == 0)
	// {
	// 	ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc",
	// 		"Unrecognized functional type in '" + xc_func_in + "'.");
	// }

    // determine the functional id
	std::vector<int> func_id(xcfunc_words_.size(), -1);
	std::transform(xcfunc_words_.begin(), xcfunc_words_.end(), func_id.begin(),
		[](const std::string& func) { return xc_functional_get_number(func.c_str()); });
	// if there is any -1, it means the functional is not recognized
	const bool not_recognized_xc = std::any_of(func_id.begin(), func_id.end(), 
		[](int id) { return id == -1; });
	if (not_recognized_xc)
	{
		std::string message = "Unrecognized exchange-correlation functional '" + xc_func_in + "'.\n"
							  " Possible source: Pseudopotential file or dft_functional parameter.\n"
							  " Please explicitly set dft_functional in INPUT,\n"
							  " or verify the functional name is supported.";
		ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc", message);
	}

	// return
    return std::make_pair(func_type, func_id);
}

const std::vector<double> in_built_xc_func_ext_params(const int id)
{
	switch(id)
	{
		// finite temperature XC functionals
		case XC_LDA_XC_KSDT:
			return {PARAM.inp.xc_temperature * 0.5};
		case XC_LDA_XC_CORRKSDT:
			return {PARAM.inp.xc_temperature * 0.5};
		case XC_LDA_XC_GDSMFB:
			return {PARAM.inp.xc_temperature * 0.5};
#ifdef __EXX
		// hybrid functionals
		case XC_HYB_GGA_XC_PBEH:
			return {GlobalC::exx_info.info_global.hybrid_alpha,
					GlobalC::exx_info.info_global.hse_omega, 
					GlobalC::exx_info.info_global.hse_omega};
		case XC_HYB_GGA_XC_HSE06:
			return {GlobalC::exx_info.info_global.hybrid_alpha,
					GlobalC::exx_info.info_global.hse_omega, 
					GlobalC::exx_info.info_global.hse_omega};
		// short-range of B88_X
		case XC_GGA_X_ITYH:
			return {GlobalC::exx_info.info_global.hse_omega};
		// short-range of LYP_C
		case XC_GGA_C_LYPR:
			return {0.04918, 0.132, 0.2533, 0.349, 
					0.35/2.29, 2.0/2.29, GlobalC::exx_info.info_global.hse_omega};
#endif
		default:
			return std::vector<double>{};
	}
}

const std::vector<double> external_xc_func_ext_params(const int id)
{
	const std::map<int, std::vector<double>> mymap = {
		{
			PARAM.inp.xc_exch_ext[0],
			std::vector<double>(PARAM.inp.xc_exch_ext.begin()+1,
								PARAM.inp.xc_exch_ext.end())
		},
		{
			PARAM.inp.xc_corr_ext[0],
			std::vector<double>(PARAM.inp.xc_corr_ext.begin()+1,
								PARAM.inp.xc_corr_ext.end())
		}
 	};
	auto it = mymap.find(id);
	return (it != mymap.end()) ? it->second : std::vector<double>{};
}

std::vector<xc_func_type> 
XC_Functional_Libxc::init_func(const std::vector<int> &func_id, 
							   const int xc_polarized)
{
	std::vector<xc_func_type> funcs;
	for (int id : func_id)
	{
		funcs.push_back({}); // create placeholder
		xc_func_init(&funcs.back(), id, xc_polarized); // instantiate the XC term

		// search for external parameters
		const std::vector<double> in_built_ext_params = in_built_xc_func_ext_params(id);
		const std::vector<double> external_ext_params = external_xc_func_ext_params(id);
		// for temporary use, I name their size as n1 and n2
		const int n1 = in_built_ext_params.size();
		const int n2 = external_ext_params.size();

// #ifdef __DEBUG // will the following assertion cause performance issue?
		// assert the number of parameters should be either zero or the value from
		// libxc function xc_func_info_get_n_ext_params, this is to avoid the undefined
		// behavior of illegal memory access
		const xc_func_info_type* info = xc_func_get_info(&funcs.back());
		const int nref = xc_func_info_get_n_ext_params(info);
		assert ((n1 == 0) || (n1 == nref) || (n2 == 0) || (n2 == nref));
// #endif

		// external overwrites in-built if the same functional id is found in both maps
		const double* xc_func_ext_params = 
			(n2 > 0) ? external_ext_params.data() : 
			(n1 > 0) ? in_built_ext_params.data() :
			nullptr; // nullptr if no external parameters are found

		// if there are no external parameters, do nothing, otherwise we set
		if(xc_func_ext_params != nullptr)
		{
			// set the external parameters
			xc_func_set_ext_params(&funcs.back(), const_cast<double*>(xc_func_ext_params));
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