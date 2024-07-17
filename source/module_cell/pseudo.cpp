#include "pseudo.h"
#include "module_base/tool_title.h"

pseudo::pseudo()
{
}

pseudo::~pseudo()
{
    delete[] els;
    delete[] lchi;
    delete[] oc;
    delete[] jjj;
    delete[] jchi;
    delete[] nn;
    delete[] vloc_at;
    delete[] r;
    delete[] rab;
    delete[] rho_at;
    delete[] rho_atc;
    delete[] lll;
}

//---------------------------------------------------------------------
void pseudo::set_pseudo()
{
    ModuleBase::TITLE("pseudo", "set_pseudo");

    // call subroutines
    this->set_pseudo_h();
	this->set_pseudo_atom();
	this->set_pseudo_vl();

	if (nbeta == 0) {
		return;
}
	
	if (lll == nullptr)
	{
		lll = new int[nbeta];
		assert(lll != nullptr);
	}

	nh = 0;

	for (int nb = 0; nb < nbeta;nb++)
	{
		nh += 2 * lll [nb] + 1;
	} 

	return;
} // end subroutine set_pseudo_upf

void pseudo::print_pseudo(std::ofstream& ofs)
{
	print_pseudo_vl(ofs);
	ofs << "\n pseudo : ";
	ofs << "\n kkbeta	" << kkbeta;
	ofs << "\n nh  " << nh;
	output::printr1_d(ofs, " lll : ", lll, nbeta);
	output::printrm(ofs, " betar : ", betar);
	output::printrm(ofs, " dion : ", dion);
	ofs << "\n ----------------------";
}


void pseudo::set_pseudo_h()
{
	ModuleBase::TITLE("pseudo","set_pseudo_h");
	
	// mohan update 2021-02-22
	//  max number of points in the atomic radial mesh
	int ndmx = 200000; 
	if (this->mesh > ndmx)
	{
		std::cout << "\n set_pseudo_h, too many grid points,";
	}

	if (this->els == nullptr)
	{
		this->els = new std::string[nchi];
		assert(this->els != nullptr);
	}

	if (this->lchi == nullptr)
	{
		this->lchi = new int[this->nchi];
		assert(this->lchi != nullptr);
	}

	if (this->oc == nullptr)
	{
		this->oc = new double[nchi];
		assert(this->oc != nullptr);
	}

	if (jjj == nullptr) {
		this->jjj = new double[nbeta];
		assert(this->jjj != nullptr);
		assert(!this->has_so);
		for (int i=0; i<nbeta; i++)
		{
			this->jjj[i]  = 0;
		}
	}

	if (nn == nullptr) {
		this->nn = new int[nchi];
		assert(this->nn != nullptr);
		assert(!this->has_so);
		for (int i=0; i<nchi; i++)
		{
			this->nn[i] = 0;
		}
	}

	if (jchi == nullptr) {
		this->jchi = new double[nchi];
		assert(this->jchi != nullptr);
		assert(!this->has_so);
		for (int i=0; i<nchi; i++)
		{
			this->jchi[i] = 0;
		}
	}

    return;
} // end subroutine set_pseudo_upf


void pseudo::set_pseudo_atom()
{
	ModuleBase::TITLE("pseudo","set_pseudo_atom");

	// mohan 2009-12-15
	// mohan update again 2011-05-23, 
	// in order to calculate more accurate Vna.
	this->rcut = GlobalV::PSEUDORCUT;//(a.u.);
	
// remember to update here if you need it.
//	rcut = 25.0; 

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"PAO radial cut off (Bohr)",rcut);
	if(rcut <= 0.0)
	{
		ModuleBase::WARNING_QUIT("pseudo_atom::set_pseudo_atom","PAO rcut<=0.0");
	}

	// chi.create(nchi, mesh);

	if (r == nullptr) {
		r = new double[mesh];
		assert(r != nullptr);
		ModuleBase::GlobalFunc::ZEROS(r, mesh);
	}

	if (rab == nullptr) {
		rab = new double[mesh];
		assert(rab != nullptr);
		ModuleBase::GlobalFunc::ZEROS(rab, mesh);
	}

	if (rho_at == nullptr) {
		rho_at = new double[mesh];
		assert(rho_at != nullptr);
		ModuleBase::GlobalFunc::ZEROS(rho_at, mesh);
	}

	if (rho_atc == nullptr) {
		rho_atc = new double[mesh];
		assert(rho_atc != nullptr);
		assert(!nlcc);
		ModuleBase::GlobalFunc::ZEROS(rho_atc, mesh);
	}

	bool br = false;

	this->msh = 0;

	for (int ir = 0;ir < mesh;ir++)
	{
		if (r [ir] > rcut)
		{
			msh = ir + 1;
			br = true;
			break;
		}
	}

	if (br)
	{
		// force msh to be odd for simpson integration
		msh = 2 * (int)((msh + 1) / 2) - 1;	// 5
	}
	else
	{
		msh = mesh ;
	}

	return;
} // end subroutine set_pseudo



void pseudo::set_pseudo_vl()
{
	ModuleBase::TITLE("pseudo","set_pseudo_vl");

	assert(mesh>0);//mohan add 2021-05-01

	if (vloc_at == nullptr) {
		vloc_at = new double[mesh];
		assert(vloc_at != nullptr);
	}


	return;
} 


void pseudo::print_pseudo_atom(std::ofstream &ofs)
{
	print_pseudo_h(ofs);
	ofs << "\n pseudo_atom : ";
	ofs << "\n msh	" << msh;
//	ofs	<< "\n nchi	" << nchi;
	output::printr1_d(ofs, " r : ", r, mesh);
	output::printr1_d(ofs, " rab : ", rab, mesh);
	output::printr1_d(ofs, " rho_atc : ", rho_atc, mesh);
	output::printr1_d(ofs, " rho_at : ", rho_at, mesh);
	output::printr1_d(ofs," jchi : ", jchi, nchi);
	output::printrm(ofs, " chi : ", chi);
	ofs << "\n ----------------------";
}


void pseudo::print_pseudo_vl(std::ofstream &ofs)
{
	ofs << "\n pseudo_vl:";
	print_pseudo_atom(ofs);
	output::printr1_d(ofs, "vloc_at : ", vloc_at, mesh);
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
    output::printr1_d(ofs, " lchi: ", lchi, nchi);
    output::printr1_d(ofs, " oc: ", oc, nchi);
    ofs << "\n ----------------------";
}

