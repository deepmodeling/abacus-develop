#include "xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

#ifdef USE_LIBXC
#include "xc_functional_libxc.h"
#endif

XC_Functional::XC_Functional(){}

XC_Functional::~XC_Functional(){}

std::vector<int> XC_Functional::func_id(1);
int XC_Functional::func_type = 0;
bool XC_Functional::use_libxc = true;
double XC_Functional::hybrid_alpha = 0.25;
std::map<int, double> XC_Functional::scaling_factor_xc = { {1, 1.0} }; // added by jghan, 2024-10-10

void XC_Functional::set_hybrid_alpha(const double alpha_in)
{
    hybrid_alpha = alpha_in;
}

double XC_Functional::get_hybrid_alpha()
{
    return hybrid_alpha;
}

int XC_Functional::get_func_type()
{
    return func_type;
}
void XC_Functional::set_xc_first_loop(const UnitCell& ucell)
{
    /** In the special "two-level" calculation case,
the first scf iteration only calculate the functional without exact
exchange. but in "nscf" calculation, there is no need of "two-level"
method. */
    if (ucell.atoms[0].ncpp.xc_func == "HF"
        || ucell.atoms[0].ncpp.xc_func == "PBE0"
        || ucell.atoms[0].ncpp.xc_func == "HSE") {
        XC_Functional::set_xc_type("pbe");
    }
    else if (ucell.atoms[0].ncpp.xc_func == "SCAN0") {
        XC_Functional::set_xc_type("scan");
    }
}

// The setting values of functional id according to the index in LIBXC
// for detail, refer to https://www.tddft.org/programs/libxc/functionals/
void XC_Functional::set_xc_type(const std::string xc_func_in)
{
    //Note : due to the separation of gcx_spin and gcc_spin,
    //when you are adding new GGA functionals,
    //please put exchange first, followed by correlation,
    //such as for PBE we have:
    //        func_id.push_back(XC_GGA_X_PBE);
    //        func_id.push_back(XC_GGA_C_PBE);

    func_id.clear();
    scaling_factor_xc.clear(); // added by jghan, 2024-07-07
    std::string xc_func = xc_func_in;
    std::transform(xc_func.begin(), xc_func.end(), xc_func.begin(), (::toupper));
	if( xc_func == "LDA" || xc_func == "PZ" || xc_func == "SLAPZNOGXNOGC") //SLA+PZ
	{
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PZ);
        func_type = 1;
        use_libxc = false;
#ifdef USE_PAW
        if(PARAM.inp.use_paw)
        {
            if(PARAM.inp.nspin != 1)
            {
                ModuleBase::WARNING_QUIT("set_xc_type","paw does not support pz with spin polarization");
            }
            else
            {
                GlobalC::paw_cell.set_libpaw_xc(1,2);
            }
        }
#endif
	}
    else if (xc_func == "PWLDA")
    {
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PW);
        func_type = 1;
        use_libxc = false;
#ifdef USE_PAW
        if(PARAM.inp.use_paw) { GlobalC::paw_cell.set_libpaw_xc(1,7);
}
#endif
    }
	else if ( xc_func == "PBE" || xc_func == "SLAPWPBXPBC") //PBX+PBC
	{
        func_id.push_back(XC_GGA_X_PBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
#ifdef USE_PAW
        if(PARAM.inp.use_paw) { GlobalC::paw_cell.set_libpaw_xc(2,11);
}
#endif
	}
	else if ( xc_func == "PBESOL") //PBX_S+PBC_S
	{
        func_id.push_back(XC_GGA_X_PBE_SOL);
        func_id.push_back(XC_GGA_C_PBE_SOL);
        func_type = 2;
        use_libxc = false;
	}
	else if( xc_func == "REVPBE" ) //PBX_r+PBC
	{
		func_id.push_back(XC_GGA_X_PBE_R);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
#ifdef USE_PAW
        if(PARAM.inp.use_paw) { GlobalC::paw_cell.set_libpaw_xc(2,14);
}
#endif
	}
	else if ( xc_func == "WC") //WC+PBC
	{
        func_id.push_back(XC_GGA_X_WC);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "BLYP") //B88+LYP
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "BP") //B88+P86
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_P86);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "PW91") //PW91_X+PW91_C
	{
        func_id.push_back(XC_GGA_X_PW91);
        func_id.push_back(XC_GGA_C_PW91);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "HCTH") //HCTH_X+HCTH_C
	{
        func_id.push_back(XC_GGA_X_HCTH_A);
        func_id.push_back(XC_GGA_C_HCTH_A);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "OLYP") //OPTX+LYP
	{
        func_id.push_back(XC_GGA_X_OPTX);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
#ifdef USE_LIBXC
	else if ( xc_func == "SCAN")
	{
        func_id.push_back(XC_MGGA_X_SCAN);
        func_id.push_back(XC_MGGA_C_SCAN);
        func_type = 3;
        use_libxc = true;
	}
    else if ( xc_func == "SCAN0")
	{
        func_id.push_back(XC_MGGA_X_SCAN);
        func_id.push_back(XC_MGGA_C_SCAN);
        func_type = 5;
        use_libxc = true;
	}
#endif
    else if( xc_func == "HF")
    {
        func_type = 4;
        use_libxc = false;
    }
   	else if( xc_func == "PBE0")
	{
        func_id.push_back(XC_HYB_GGA_XC_PBEH);
        func_type = 4;
        use_libxc = false;
	}
    else if( xc_func == "OPT_ORB" ||  xc_func == "NONE" || xc_func == "NOX+NOC")
    {
        // not doing anything
    }
    else if( xc_func == "MULLER" || xc_func == "POWER" ) // added by jghan, 2024-07-06
    {
        func_type = 4;
        use_libxc = false;
    }
#ifdef USE_LIBXC
    else if( xc_func == "HSE")
    {
        func_id.push_back(XC_HYB_GGA_XC_HSE06);
        func_type = 4;
        use_libxc = true;
    }
    // added by jghan, 2024-07-06
    else if( xc_func == "WP22")
    {
        func_id.push_back(XC_GGA_X_ITYH);   // short-range of B88_X, id=529
        func_id.push_back(XC_GGA_C_LYPR);   // short-range of LYP_C, id=624
        func_type = 4;
        use_libxc = true;
    }
    else if( xc_func == "CWP22")
    {   
        // BLYP_XC_lr = -BLYP_XC_sr + BLYP_XC, the realization of it is in v_xc_libxc() function, xc_functional_libxc_vxc.cpp
        func_id.push_back(XC_GGA_X_ITYH);   // short-range of B88_X, id=529
        func_id.push_back(XC_GGA_C_LYPR);   // short-range of LYP_C, id=624
        func_id.push_back(XC_GGA_X_B88);    // complete B88_X, id=106
        func_id.push_back(XC_GGA_C_LYP);    // complete LYP_C, id=131

        // the scaling factor of CWP22-functionals
        scaling_factor_xc[XC_GGA_X_ITYH] = -1.0;
        scaling_factor_xc[XC_GGA_C_LYPR] = -1.0;
        scaling_factor_xc[XC_GGA_X_B88] = 1.0;
        scaling_factor_xc[XC_GGA_X_B88] = 1.0;

        func_type = 4;
        use_libxc = true;
    }
    else if( xc_func == "BLYP_LR")
    {   
        // BLYP_XC_lr = -BLYP_XC_sr + BLYP_XC, the realization of it is in v_xc_libxc() function, xc_functional_libxc_vxc.cpp
        func_id.push_back(XC_GGA_X_ITYH);   // short-range of B88_X, id=529
        func_id.push_back(XC_GGA_C_LYPR);   // short-range of LYP_C, id=624
        func_id.push_back(XC_GGA_X_B88);    // complete B88_X, id=106
        func_id.push_back(XC_GGA_C_LYP);    // complete LYP_C, id=131

        // the scaling factor of BLYP_LR-functionals
        scaling_factor_xc[XC_GGA_X_ITYH] = -1.0;
        scaling_factor_xc[XC_GGA_C_LYPR] = -1.0;
        scaling_factor_xc[XC_GGA_X_B88] = 1.0;
        scaling_factor_xc[XC_GGA_X_B88] = 1.0;

        func_type = 2;
        use_libxc = true;
    }
#endif
    else
    {
#ifdef USE_LIBXC
        //see if it matches libxc functionals
        const std::pair<int,std::vector<int>> type_id = XC_Functional_Libxc::set_xc_type_libxc(xc_func);
        func_type = std::get<0>(type_id);
        func_id = std::get<1>(type_id);
        use_libxc = true;
#else
        ModuleBase::WARNING_QUIT("xc_functional.cpp","functional name not recognized!");
#endif
    }

	if (func_id[0] == XC_GGA_X_OPTX)
	{
		std::cerr << "\n OPTX untested please test,";
	}

    if((func_type == 4 || func_type == 5) && PARAM.inp.basis_type == "pw")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","hybrid functional not realized for planewave yet");
    }
    if((func_type == 3 || func_type == 5) && PARAM.inp.nspin==4)
    {
        ModuleBase::WARNING_QUIT("set_xc_type","meta-GGA has not been implemented for nspin = 4 yet");
    }
    //if((func_type == 3 || func_type == 5) && PARAM.inp.cal_stress == 1 && PARAM.inp.nspin!=1)
    //{
    //    ModuleBase::WARNING_QUIT("set_xc_type","mgga stress not implemented for polarized case yet");
    //}

#ifndef __EXX
    if(func_type == 4 || func_type == 5)
    {
        ModuleBase::WARNING_QUIT("set_xc_type","compile with libri to use hybrid functional");
    }
#endif

#ifndef USE_LIBXC
    if(xc_func == "SCAN" || xc_func == "HSE" || xc_func == "SCAN0" 
        || xc_func == "MULLER" || xc_func == "POWER" || xc_func == "WP22" || xc_func == "CWP22")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","to use SCAN, SCAN0, or HSE, LIBXC is required");
    }
    use_libxc = false;
#endif

}
