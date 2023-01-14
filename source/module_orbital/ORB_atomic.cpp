#include "ORB_atomic.h"

Numerical_Orbital_AtomRelation Numerical_Orbital::NOAR;

Numerical_Orbital::Numerical_Orbital()
{
	//make std::pair of new and delete
	//question remains
	this->nchi = nullptr;
	this->phiLN = new Numerical_Orbital_Lm[1];
}

Numerical_Orbital::~Numerical_Orbital()
{
	delete[] nchi;
	delete[] phiLN;
}

void Numerical_Orbital::set_orbital_info
(
    const int& type_in,
    const std::string& label_in,
    const int& lmax_in,
    const int* nchi_in,
    const int& total_nchi_in
)
{
	//what is GlobalV::test_overlap
	ModuleBase::TITLE("Numerical_Orbital", "set_type_info");

	// (1) set type,label,lmax
	this->type = type_in;
	this->label = label_in;
	this->lmax = lmax_in;

	// (2) set nchi and total nchi.
	delete[] this->nchi;
	this->nchi = new int[this->lmax+1];
	for (int i = 0; i < this->lmax + 1; i++)
	{
		this->nchi[i] = nchi_in[i];
	}

	// we need this to generate numerical_orbital_lm.
	if (total_nchi_in < 0 || total_nchi_in > 50) 
	{
		ModuleBase::WARNING_QUIT("Numerical_Orbital::init", "total_nchi<0 || total_nchi>50");
	}
	else 
	{
		this->total_nchi = total_nchi_in;
	}

	// (3) set the rcut and check the rcut
	this->rcut = 0.0;
	for (int i=0; i< total_nchi_in; i++)
	{
		this->rcut = this->phiLN[i].rcut;
		for(int j=0; j< total_nchi_in; j++)
		{
			assert( rcut == this->phiLN[j].rcut );
		}
	}
	assert(rcut > 0.0);

	// (4) set max_nchi
	this->max_nchi=0;
	for(int L=0; L<lmax+1; L++)
	{
		max_nchi = std::max( max_nchi, nchi[L] ); 
	}

	// (8) set find_chi
	assert( lmax+1 > 0 );
	this->find_chi.create( lmax+1, max_nchi );
	int ichi=0;
	for(int L=0; L<=lmax; ++L)
	{
		for(int N=0; N<nchi[L]; ++N)
		{
			find_chi(L,N) = ichi;
			++ichi;
		}
	}
	assert(ichi == total_nchi);

	return;
}
