#include "pseudo.h"
#include "source_base/tool_title.h"

pseudo::pseudo()
{
}

pseudo::~pseudo()
{
}

void pseudo::check_betar()
{
	bool min_flag = false;
	for (int ib = 0; ib < nbeta; ib++)
	{
		for (int ir = 0; ir < mesh; ir++)
		{
			// Get the bit representation of the double
			uint64_t bits = *(uint64_t*)&betar(ib, ir);
    		// Extract exponent field (bits 52-62)
			uint64_t exponent = (bits >> 52) & 0x7FF;
			// Define exponent threshold for 1e-30
    		// Calculated as: bias + floor(log2(1e-30))
    		// Where bias = 1023 and log2(1e-30) ≈ -99.657
			// Thus threshold is approximately 923
			if ((exponent <= 923))
			{
				min_flag = true;
				betar(ib, ir) = 0.0;
			}
		}
	}
	if (min_flag)
	{
		std::cout << "WARNING: some of potential function is set to zero cause of less than 1e-30.\n";
	}
}

void pseudo::print_pseudo(std::ofstream& ofs)
{
	print_pseudo_vl(ofs);
	ofs << "\n pseudo : ";
	ofs << "\n kkbeta	" << kkbeta;
	ofs << "\n nh  " << nh;
	output::printr1_d(ofs, " lll : ", lll.data(), nbeta);
	output::printrm(ofs, " betar : ", betar);
	output::printrm(ofs, " dion : ", dion);
	ofs << "\n ----------------------";
}

void pseudo::print_pseudo_atom(std::ofstream &ofs)
{
	print_pseudo_h(ofs);
	ofs << "\n pseudo_atom : ";
	ofs << "\n msh	" << msh;
//	ofs	<< "\n nchi	" << nchi;
	output::printr1_d(ofs, " r : ", r.data(), mesh);
	output::printr1_d(ofs, " rab : ", rab.data(), mesh);
	output::printr1_d(ofs, " rho_atc : ", rho_atc.data(), mesh);
	output::printr1_d(ofs, " rho_at : ", rho_at.data(), mesh);
	output::printr1_d(ofs," jchi : ", jchi.data(), nchi);
	output::printrm(ofs, " chi : ", chi);
	ofs << "\n ----------------------";
}


void pseudo::print_pseudo_vl(std::ofstream &ofs)
{
	ofs << "\n pseudo_vl:";
	print_pseudo_atom(ofs);
	output::printr1_d(ofs, "vloc_at : ", vloc_at.data(), mesh);
	ofs << "\n ----------------------------------- ";
}

void pseudo::print_pseudo_h(std::ofstream &ofs)
{
    ofs << "\n pseudo_info :";
    ofs << "\n nv       " << nv;
    ofs << "\n psd  " << psd;
    ofs << "\n pp_type  " << pp_type;
    ofs << "\n tvanp    " << tvanp;
    ofs << "\n nlcc " << nlcc;
    ofs << "\n dft  " << xc_func;
    ofs << "\n zv       " << zv;
    ofs << "\n etotps   " << etotps;
    ofs << "\n ecutwfc " << ecutwfc;
    ofs << "\n ecutrho " << ecutrho;
    ofs << "\n lmax " << lmax;
    ofs << "\n mesh " << mesh;
    ofs << "\n nchi " << nchi;
    ofs << "\n nbeta    " << nbeta;
//  out.printr1_d(ofs," els: ", els, nchi);
    output::printr1_d(ofs, " lchi: ", lchi.data(), nchi);
    output::printr1_d(ofs, " oc: ", oc.data(), nchi);
    ofs << "\n ----------------------";
}

