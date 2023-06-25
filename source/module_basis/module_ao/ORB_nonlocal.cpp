#include "ORB_nonlocal.h"
#include "module_base/global_function.h"

Numerical_Nonlocal::Numerical_Nonlocal()
{
	//make std::pair of new and delete
	//question remains
	this->type = 0;
	this->lmax = 0;
	this->Proj = new Numerical_Nonlocal_Lm[1];
	this->nproj = -1;
	//zhengdy-soc, for optimize nonlocal part
}

Numerical_Nonlocal::~Numerical_Nonlocal()
{
	delete[] Proj;
}

void Numerical_Nonlocal::set_type_info
(
	const int& type_in,
	const std::string& label_in,
	const std::string& type_ps_in,
	const int& lmax_in,
	const int& nproj_in,
	const Numerical_Nonlocal_Lm* Proj_in
)
{
	//ModuleBase::TITLE("Numerical_Nonlocal","set_type_info");

	this->type = type_in;
	this->label = label_in;
	this->type_ps = type_ps_in;

	if (lmax_in < -1 || lmax_in > 20)
	{
		ModuleBase::WARNING_QUIT("Numerical_Nonlocal", "bad input of lmax : should be between -1 and 20");
	}

	this->lmax = lmax_in;

	this->nproj = nproj_in;

	assert(nproj >= 0);

//----------------------------------------------------------
// EXPLAIN : non_local pseudopotential projector for each l
//----------------------------------------------------------
	//only store radial function
	delete[] Proj;
	this->Proj = new Numerical_Nonlocal_Lm[this->nproj];

	for (int p1=0; p1<nproj; p1++)
	{
		this->Proj[p1] = Proj_in[p1];
	}

	this->rcut_max = 0.0;
	for(int p1=0; p1<nproj; p1++)
	{
		this->rcut_max = std::max( this->Proj[p1].getRcut(), rcut_max ); 
	}
	return;
}

