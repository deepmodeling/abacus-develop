#include <cassert>
#include <iostream>
#include <sstream>
#include <cmath>
#include "ORB_nonlocal_lm.h"
#include "../src_global/math_integral.h"
#include "../src_global/global_function.h"
#include "../src_global/mathzone.h" // use Polynomial_Interpolation_xy
#include "../src_global/mathzone_add1.h" // use SplineD2

Numerical_Nonlocal_Lm::Numerical_Nonlocal_Lm()
{
	label = "";
	index_atom_type = 0;
	angular_momentum_l = 0;
	
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

void Numerical_Nonlocal_Lm::renew(void)
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
	ZEROS(r_radial, nr);
	ZEROS(rab, nr);
	ZEROS(beta_r, nr);
	ZEROS(beta_uniform, nr_uniform);
	ZEROS(dbeta_uniform, nr_uniform);
	ZEROS(k_radial, nk);
	ZEROS(beta_k, nk);
}

void Numerical_Nonlocal_Lm::freemem(void)
{
	delete[] this->r_radial;
	delete[] this->rab;
	delete[] this->beta_r;
	delete[] this->beta_uniform;
	delete[] this->dbeta_uniform;
	delete[] this->k_radial;
	delete[] this->beta_k;
}

Numerical_Nonlocal_Lm& Numerical_Nonlocal_Lm::operator=
(
    const Numerical_Nonlocal_Lm & nol
)
{
	this->label = nol.label;
	this->index_atom_type = nol.index_atom_type;
	this->angular_momentum_l = nol.angular_momentum_l;

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
 	const string &label_in,
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

void Numerical_Nonlocal_Lm::extra_uniform(const double &dr_uniform_in)
{
	assert(dr_uniform_in>0.0);
	this->dr_uniform = dr_uniform_in;
	this->nr_uniform = static_cast<int>(rcut/dr_uniform) + 10;

//	cout << " nr_uniform = " << nr_uniform << endl;
	
	delete[] this->beta_uniform;
	this->beta_uniform = new double[nr_uniform];
	ZEROS (this->beta_uniform, nr_uniform);

	// mohan fix bug 2011-04-14
	// the beta_r is beta*r.
	// and beta is beta!!!
	double* beta = new double[nr];
	ZEROS(beta,nr);
	for(int ir=0; ir<nr; ir++)
	{
		assert(r_radial[ir]>0.0);
		beta[ir] = this->beta_r[ir]/r_radial[ir];
	}

	for (int ir = 0; ir < this->nr_uniform; ir++)
	{
		double rnew = ir * dr_uniform;
		this->beta_uniform[ir] = Mathzone::Polynomial_Interpolation_xy(this->r_radial, beta, this->nr, rnew); 
    }
	delete[] beta;

	delete[] this->dbeta_uniform;
	this->dbeta_uniform = new double[nr_uniform];
	ZEROS(this->dbeta_uniform, nr_uniform);
	
	double* y2 = new double[nr];
	Mathzone_Add1::SplineD2 (r_radial, beta_r, nr, 0.0, 0.0, y2);

	double* rad = new double[nr_uniform];
	for (int ir = 0; ir < nr_uniform; ir++)
	{
		rad[ir] = ir*dr_uniform;
	}
	
	double* tmp = new double[nr_uniform];
	Mathzone_Add1::Cubic_Spline_Interpolation(r_radial, beta_r, y2, 
					nr, rad, nr_uniform, tmp, dbeta_uniform );

	for(int ir= 0 ; ir<nr_uniform; ir++)
	{
		//assert(dbeta_uniform[ir]==beta_uniform[ir]);
	}

	delete [] y2;
	delete [] rad;	
	delete [] tmp;
	return;
}


void Numerical_Nonlocal_Lm::get_kradial(void)
{
    //TITLE("Numerical_Nonlocal_Lm","get_kradial");
    double *jl = new double[nr];
    double *integrated_func = new double[nr];

	// mohan add constant of pi, 2021-04-26
	const double pi = 3.14159265358979323846;
    const double pref = sqrt( 2.0 / pi );

    for (int ik = 0; ik < nk; ik++)
    {
        Mathzone::Spherical_Bessel(
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

        Integral::Simpson_Integral(
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
	string orbital_type;
	switch( this->angular_momentum_l )
	{
		case 0: orbital_type = "s"; break;
		case 1: orbital_type = "p"; break;
		case 2: orbital_type = "d"; break;
		case 3: orbital_type = "f"; break;
		case 4: orbital_type = "g"; break;
		default: WARNING_QUIT("Numerical_Orbital_Lm::plot","Please check in functoin.");
	}

	if(my_rank==0)
	{
		stringstream ssr, ssk, ssru;
		ssr << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-r.dat";

		ssk << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-k.dat";

		ssru << global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << "-proj-ru.dat";

		ofstream ofsr(ssr.str().c_str());
		ofstream ofsk(ssk.str().c_str());
		ofstream ofsru(ssru.str().c_str());

		if (!ofsk || !ofsr || !ofsru)
		{
			WARNING_QUIT("Numerical_Orbital_Lm::plot", "Can't open files!");
		}

		for (int i = 0; i < this->nr; i++)
		{
			ofsr << this->r_radial[i] << " " << this->beta_r[i] << endl;
		}

		for (int i = 0; i < this->nk; i++)
		{
			ofsk << this->k_radial[i] << " " << this->beta_k[i] << endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsru << i * this->dr_uniform << " " << this->beta_uniform[i] << endl;
		}
		
		ofsr.close();
		ofsk.close();
		ofsru.close();
	}
	return;
}
