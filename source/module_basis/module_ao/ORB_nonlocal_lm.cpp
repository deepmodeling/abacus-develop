#include "ORB_nonlocal_lm.h"

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/math_integral.h"
#include "module_base/math_polyint.h"
#include "module_base/math_sphbes.h"
#include "module_base/mathzone.h"      /// use Polynomial_Interpolation_xy, Spherical_Bessel
#include "module_base/mathzone_add1.h" /// use SplineD2
#include "module_parameter/parameter.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

Numerical_Nonlocal_Lm::Numerical_Nonlocal_Lm()
{
	label = "";
	index_atom_type = 0;
	angular_momentum_l = 0;
	index_proj = 0;
	
	nr = 1;
	nk = 1;

	rcut = 0.0;
	kcut = 0.0;
	dk = 0.0;

	nr_uniform = 1;
	dr_uniform = -1.0;
	this->renew();
}

Numerical_Nonlocal_Lm::~Numerical_Nonlocal_Lm()
{
	this->freemem();
}

void Numerical_Nonlocal_Lm::renew()
{
	assert(nr_uniform>0);
	assert(nr>0);
	assert(nk>0);
	this->r_radial = new double[nr];
	this->rab = new double[nr];
	this->beta_r = new double[nr];
	this->beta_uniform = new double[nr_uniform];
	this->dbeta_uniform = new double[nr_uniform];
	this->k_radial = new double[nk];
	this->beta_k = new double[nk];
	ModuleBase::GlobalFunc::ZEROS(r_radial, nr);
	ModuleBase::GlobalFunc::ZEROS(rab, nr);
	ModuleBase::GlobalFunc::ZEROS(beta_r, nr);
	ModuleBase::GlobalFunc::ZEROS(beta_uniform, nr_uniform);
	ModuleBase::GlobalFunc::ZEROS(dbeta_uniform, nr_uniform);
	ModuleBase::GlobalFunc::ZEROS(k_radial, nk);
	ModuleBase::GlobalFunc::ZEROS(beta_k, nk);
}

void Numerical_Nonlocal_Lm::freemem()
{
	delete[] this->r_radial;
	delete[] this->rab;
	delete[] this->beta_r;
	delete[] this->beta_uniform;
	delete[] this->dbeta_uniform;
	delete[] this->k_radial;
	delete[] this->beta_k;

    r_radial = nullptr;
    rab = nullptr;
    beta_r = nullptr;
    beta_uniform = nullptr;
    dbeta_uniform = nullptr;
    k_radial = nullptr;
    beta_k = nullptr;
}

Numerical_Nonlocal_Lm& Numerical_Nonlocal_Lm::operator=
(
    const Numerical_Nonlocal_Lm & nol
)
{
	this->label = nol.label;
	this->index_atom_type = nol.index_atom_type;
	this->angular_momentum_l = nol.angular_momentum_l;
    this->index_proj = nol.index_proj;

	this->nr = nol.nr;
	this->nk = nol.nk;

	this->nr_uniform = nol.nr_uniform;
	this->dr_uniform = nol.dr_uniform;

	this->rcut = nol.rcut;
	this->kcut = nol.kcut;

	this->dk = nol.dk;

	this->freemem();
	this->renew();

	for (int ir = 0; ir < nol.nr; ir++)
	{
		this->r_radial[ir] = nol.r_radial[ir];
		this->rab[ir] = nol.rab[ir];
		this->beta_r[ir] = nol.beta_r[ir];
	}

	for (int ir = 0; ir < nr_uniform; ir++)
	{
		this->beta_uniform[ir] = nol.beta_uniform[ir];
		this->dbeta_uniform[ir] = nol.dbeta_uniform[ir];
	}

	for (int ik = 0; ik < nol.nk; ik++)
	{
		this->k_radial[ik] = nol.k_radial[ik];
		this->beta_k[ik] = nol.beta_k[ik];
	}

	return *this;
}

void Numerical_Nonlocal_Lm::set_NL_proj(
 	const std::string &label_in,
    const int &index_atom_type_in,
    const int &angular_momentum_l_in,
    const int &nr_in,
    const double *rab_in,
    const double *r_radial_in,
    const double *beta_r_in,
    const int &nk_in,
    const double &dk_in,
	const double &dr_uniform_in)
{
	this->label = label_in;
	this->index_atom_type = index_atom_type_in;
	
	this->angular_momentum_l = angular_momentum_l_in;
	assert(angular_momentum_l_in>=-1); // -1 means no angular momentum.

	this->dr_uniform = dr_uniform_in;
	
	this->nr = nr_in;
	assert(nr_in>1 && nr_in <10000);
	assert(nr%2!=0);
    
	this->rcut = r_radial_in[nr-1];
    assert(rcut>=0.0);

	this->nk = nk_in;
	assert(nk%2!=0);
    
	this->dk = dk_in;
	assert(dk>0.0);

    this->freemem();
    this->renew();
	
	for (int ir = 0; ir < nr; ir++)
	{
		this->r_radial[ir] = r_radial_in[ir];
		this->rab[ir] = rab_in[ir];
		this->beta_r[ir] = beta_r_in[ir];
	}

	for (int ik = 0; ik < nk; ik++)
	{
		this->k_radial[ik] = ik * this->dk;
	}
	this->kcut = (nk-1) * this->dk;

	// (1) extra the uniform mesh
	//this->extra_uniform(dr_uniform);
	// (2) get the beta_k
	this->get_kradial();

	return;
}

void Numerical_Nonlocal_Lm::get_kradial()
{
    //ModuleBase::TITLE("Numerical_Nonlocal_Lm","get_kradial");
    double *jl = new double[nr];
    double *integrated_func = new double[nr];

    const double pref = sqrt(2.0 / ModuleBase::PI);

    for (int ik = 0; ik < nk; ik++)
    {
        ModuleBase::Sphbes::Spherical_Bessel(
                this->nr,
                this->r_radial,
                this->k_radial[ik],
                this->angular_momentum_l,
                jl);

        for (int ir = 0; ir < nr; ir++)
        {
			// beta_r is beta*r;
            integrated_func[ir] = this->beta_r[ir] * this->r_radial[ir] * jl[ir];
        }

        ModuleBase::Integral::Simpson_Integral(
                this->nr,
                integrated_func,
                this->rab,
                this->beta_k[ik]);

        this->beta_k[ik] *= ( pref*k_radial[ik]);
    }

    delete[] integrated_func;
    delete[] jl;
}


void Numerical_Nonlocal_Lm::plot(const int &my_rank)const
{
	std::string orbital_type;
	switch( this->angular_momentum_l )
	{
		case 0: orbital_type = "s"; break;
		case 1: orbital_type = "p"; break;
		case 2: orbital_type = "d"; break;
		case 3: orbital_type = "f"; break;
		case 4: orbital_type = "g"; break;
		default: ModuleBase::WARNING_QUIT("Numerical_Orbital_Lm::plot","Please check in functoin.");
	}

#ifdef __NORMAL

#else
	if(my_rank==0)
	{
		std::stringstream ssr, ssk, ssru;
		ssr << PARAM.globalv.global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-r.dat";

		ssk << PARAM.globalv.global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-k.dat";

		ssru << PARAM.globalv.global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-ru.dat";

		std::ofstream ofsr(ssr.str().c_str());
		std::ofstream ofsk(ssk.str().c_str());
		std::ofstream ofsru(ssru.str().c_str());

		if (!ofsk || !ofsr || !ofsru)
		{
			ModuleBase::WARNING_QUIT("Numerical_Orbital_Lm::plot", "Can't open files!");
		}

		for (int i = 0; i < this->nr; i++)
		{
			ofsr << this->r_radial[i] << " " << this->beta_r[i] << std::endl;
		}

		for (int i = 0; i < this->nk; i++)
		{
			ofsk << this->k_radial[i] << " " << this->beta_k[i] << std::endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsru << i * this->dr_uniform << " " << this->beta_uniform[i] << std::endl;
		}
		
		ofsr.close();
		ofsk.close();
		ofsru.close();
	}
#endif
	return;
}
