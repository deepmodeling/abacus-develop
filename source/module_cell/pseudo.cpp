#include "pseudo.h"
#include "module_base/tool_title.h"

pseudo::pseudo()
{
}

pseudo::~pseudo()
{
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

