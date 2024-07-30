#include "read_pp.h"


void Pseudopot_upf::complete_default(Atom_pseudo& pp)
{
    ModuleBase::TITLE("Pseudopot_upf", "complete_default");

    // call subroutines
    this->complete_default_h(pp);
	this->complete_default_atom(pp);
	this->complete_default_vl(pp);

	if (pp.nbeta == 0) {
		return;
	}
	
	if (pp.lll == nullptr)
	{
		pp.lll = new int[pp.nbeta];
		assert(pp.lll != nullptr);
	}

	pp.nh = 0;

	for (int nb = 0; nb < pp.nbeta;nb++)
	{
		pp.nh += 2 * pp.lll [nb] + 1;
	} 

	return;
}

void Pseudopot_upf::complete_default_h(Atom_pseudo& pp)
{
	ModuleBase::TITLE("Pseudopot_upf","complete_default_h");
	
	// mohan update 2021-02-22
	//  max number of points in the atomic radial mesh
	int ndmx = 200000; 
	if (pp.mesh > ndmx)
	{
		std::cout << "\n complete_default_h, too many grid points,";
	}

	if (pp.els == nullptr)
	{
		pp.els = new std::string[pp.nchi];
		assert(pp.els != nullptr);
	}

	if (pp.lchi == nullptr)
	{
		pp.lchi = new int[pp.nchi];
		assert(pp.lchi != nullptr);
	}

	if (pp.oc == nullptr)
	{
		pp.oc = new double[pp.nchi];
		assert(pp.oc != nullptr);
	}

	if (pp.jjj == nullptr) {
		pp.jjj = new double[pp.nbeta];
		assert(pp.jjj != nullptr);
		assert(!pp.has_so);
		for (int i=0; i<pp.nbeta; i++)
		{
			pp.jjj[i]  = 0;
		}
	}

	if (pp.nn == nullptr) {
		pp.nn = new int[pp.nchi];
		assert(pp.nn != nullptr);
		assert(!pp.has_so);
		for (int i=0; i<pp.nchi; i++)
		{
			pp.nn[i] = 0;
		}
	}

	if (pp.jchi == nullptr) {
		pp.jchi = new double[pp.nchi];
		assert(pp.jchi != nullptr);
		assert(!pp.has_so);
		for (int i=0; i<pp.nchi; i++)
		{
			pp.jchi[i] = 0;
		}
	}

    return;
}

void Pseudopot_upf::complete_default_atom(Atom_pseudo& pp)
{
	ModuleBase::TITLE("Pseudopot_upf","complete_default_atom");

	// mohan 2009-12-15
	// mohan update again 2011-05-23, 
	// in order to calculate more accurate Vna.
	pp.rcut = GlobalV::PSEUDORCUT;//(a.u.);
	
	// remember to update here if you need it.
	//	rcut = 25.0; 

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"PAO radial cut off (Bohr)",rcut);
	if(pp.rcut <= 0.0)
	{
		ModuleBase::WARNING_QUIT("Pseudopot_upf::complete_default_atom","PAO rcut<=0.0");
	}

	// chi.create(nchi, mesh);

	if (pp.r == nullptr) {
		pp.r = new double[pp.mesh];
		assert(pp.r != nullptr);
		ModuleBase::GlobalFunc::ZEROS(pp.r, pp.mesh);
	}

	if (pp.rab == nullptr) {
		pp.rab = new double[pp.mesh];
		assert(pp.rab != nullptr);
		ModuleBase::GlobalFunc::ZEROS(pp.rab, pp.mesh);
	}

	if (pp.rho_at == nullptr) {
		pp.rho_at = new double[pp.mesh];
		assert(pp.rho_at != nullptr);
		ModuleBase::GlobalFunc::ZEROS(pp.rho_at, pp.mesh);
	}

	if (pp.rho_atc == nullptr) {
		pp.rho_atc = new double[pp.mesh];
		assert(pp.rho_atc != nullptr);
		assert(!pp.nlcc);
		ModuleBase::GlobalFunc::ZEROS(pp.rho_atc, pp.mesh);
	}

	bool br = false;

	pp.msh = 0;

	for (int ir = 0;ir < pp.mesh;ir++)
	{
		if (pp.r [ir] > pp.rcut)
		{
			pp.msh = ir + 1;
			br = true;
			break;
		}
	}

	if (br)
	{
		// force msh to be odd for simpson integration
		pp.msh = 2 * (int)((pp.msh + 1) / 2) - 1;	// 5
	}
	else
	{
		pp.msh = pp.mesh ;
	}

	return;
}

void Pseudopot_upf::complete_default_vl(Atom_pseudo& pp)
{
	ModuleBase::TITLE("Pseudopot_upf","complete_default_vl");

	assert(pp.mesh>0);//mohan add 2021-05-01

	if (pp.vloc_at == nullptr) {
		pp.vloc_at = new double[pp.mesh];
		assert(pp.vloc_at != nullptr);
	}

	return;
} 