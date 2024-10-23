#include "read_pp.h"
#include "module_parameter/parameter.h"

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
	
	if (pp.lll.empty())
	{
		pp.lll = std::vector<int>(pp.nbeta, 0);
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

	if (pp.els.empty())
	{
		pp.els = std::vector<std::string>(pp.nchi, "");
	}

	if (pp.lchi.empty())
	{
		pp.lchi = std::vector<int>(pp.nchi, 0);
	}

	if (pp.oc.empty())
	{
		pp.oc = std::vector<double>(pp.nchi, 0.0);
	}

	if (pp.jjj.empty()) {
		pp.jjj = std::vector<double>(pp.nbeta, 0.0);
		assert(!pp.has_so or pp.nbeta == 0);
		for (int i=0; i<pp.nbeta; i++)
		{
			pp.jjj[i]  = 0;
		}
	}

	if (pp.nn.empty()) {
		pp.nn = std::vector<int>(pp.nchi, 0);
		assert(!pp.has_so or pp.nchi == 0);
		for (int i=0; i<pp.nchi; i++)
		{
			pp.nn[i] = 0;
		}
	}

	if (pp.jchi.empty()) {
		pp.jchi = std::vector<double>(pp.nchi, 0.0);
		assert(!pp.has_so or pp.nchi == 0);
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
	pp.rcut = PARAM.inp.pseudo_rcut;//(a.u.);
	
	// remember to update here if you need it.
	//	rcut = 25.0; 

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"PAO radial cut off (Bohr)", pp.rcut);
	if(pp.rcut <= 0.0)
	{
		ModuleBase::WARNING_QUIT("Pseudopot_upf::complete_default_atom","PAO rcut<=0.0");
	}

	// chi.create(nchi, mesh);

	if (pp.r.empty()) {
		pp.r = std::vector<double>(pp.mesh, 0.0);
	}

	if (pp.rab.empty()) {
		pp.rab = std::vector<double>(pp.mesh, 0.0);
	}

	if (pp.rho_at.empty()) {
		pp.rho_at = std::vector<double>(pp.mesh, 0.0);
	}

	if (pp.rho_atc.empty()) {
		pp.rho_atc = std::vector<double>(pp.mesh, 0.0);
		assert(!pp.nlcc or pp.mesh == 0);
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

	if (pp.vloc_at.empty()) {
		pp.vloc_at = std::vector<double>(pp.mesh, 0.0);
	}

	return;
} 