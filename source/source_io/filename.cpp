#include <set>
#include "filename.h"
#include "source_base/tool_quit.h"

namespace ModuleIO
{

std::string filename_output(
			const std::string &directory,
			const std::string &property,
			const std::string &basis,
			const int ik_local, // the ik index within each pool
			const std::vector<int> &ik2iktot,
			const int nspin,
			const int nkstot,
			const int out_type,
			const bool out_app_flag,
			const bool gamma_only,
			const int istep,
            const int iter)
{
    // output filename = "{PARAM.globalv.global_out_dir}/property{s}{spin index}
    // {k(optional)}{k-point index}{g(optional)}{geometry index1}{_basis(nao|pw)} 
    // + {".txt"/".dat"}"

	std::set<std::string> valid_properties = {"wf", "chg", "hk", "sk", "tk", "vxc"};
	if (valid_properties.find(property) == valid_properties.end()) 
	{
		ModuleBase::WARNING_QUIT("ModuleIO::filename_output", "unknown property in filename function");
	}

	std::set<std::string> valid_basis = {"pw", "nao"};
	if (valid_basis.find(basis) == valid_basis.end()) 
	{
		ModuleBase::WARNING_QUIT("ModuleIO::filename_output", "unknown basis in filename function");
	}

    assert(ik_local>=0);
    // mohan update 2025.05.07, if KPAR>1, "<" works
	assert(ik2iktot.size() <= nkstot);
    assert(nspin>0);

	// spin index
	int is0 = -1;
	// ik0 is the k-point index, starting from 0
	int ik0 = ik2iktot[ik_local];

	// spin part
	std::string spin_block;
	if(nspin == 1)
	{
	}
	else if(nspin == 2)
	{
		const int half_k = nkstot/2;
		if(ik0 >= half_k)
		{
			is0 = 2;
			ik0 -= half_k;
		}
		else
		{
			is0 = 1;
		}
        spin_block = "s" + std::to_string(is0);
	}
	else if(nspin==4)
	{
		is0 = 4;
        spin_block = "s" + std::to_string(is0);
	}


    // k-point part
    std::string kpoint_block;
    if(gamma_only)
    {
        // do nothing;
    }
    else
    {
        // mohan add 20250921
		if(property=="hk" || property=="sk" || property=="tk")
		{
			kpoint_block = std::to_string(ik0+1);
		}
		else
		{
			kpoint_block = "k" + std::to_string(ik0+1);
		}
    }

    std::string istep_block
        = (istep >= 0 && (!out_app_flag))
              ? "g" + std::to_string(istep + 1) // istep starts from 0
              : ""; // only when istep >= 0 and out_app_flag is false will write each wfc to a separate file

    std::string iter_block
        = (iter >= 0 && (!out_app_flag))
              ? "e" + std::to_string(iter) // iter starts from 1
              : ""; // only when istep >= 0 and out_app_flag is false will write each wfc to a separate file

    std::string suffix_block;
    if (out_type == 1)
    {
        suffix_block = ".txt";
    }
    else if (out_type == 2)
    {
        suffix_block = ".dat";
    }
    else
    {
        std::cout << "WARNING: the type of output wave function is not 1 or 2, so 1 is chosen." << std::endl;
        suffix_block = ".txt";
    }

    std::string fn_out
        = directory + property + kpoint_block + spin_block
          + istep_block + iter_block + "_" + basis + suffix_block;

    return fn_out;
}

}
