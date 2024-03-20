#include "Psi_Second.h"
#include "Calculate_C4.h"
#include "tools.h"

Psi_Second::Psi_Second()
{
	test=0;
	r = new double[1];
	rab = new double[1];
	psi = new double[1];
	psi_first = new double[1];
	psi_second = new double[1];
	eigen1 = new double[1];
	os_position = new double[100];
	below_max_os = new bool[1];
	count = 0;
}

Psi_Second::~Psi_Second()
{
	delete[] r;
	delete[] rab;
	delete[] psi;
	delete[] psi_first;
	delete[] psi_second;
	delete[] eigen1;
	delete[] os_position;
	delete[] below_max_os;
}

void Psi_Second::init(
	ifstream &ifs, 
	const int &ne_in, 
	const int &lmax_in,
	const double &rcut_in, 
	const int &total_nchi_in, 
	const double &tolerence)
{
	TITLE("Psi_Second","init");

	this->dr = KINETIC_DR;
	assert( dr > 0.0);
	cout << "\n KINETIC_DR=" << dr;

//	this->start = 5; // start away from 0, which may be some strange oscillation.
	this->start = 1; // mohan fix 2010-04-16
	
	// number of Jlq used.
	this->ne = ne_in;
	// lmax needed.
	this->lmax = lmax_in;
	// Cutoff radius.
	this->rcut = rcut_in;
	// Total number of radial wave functions.
	this->total_nchi = total_nchi_in;

	this->mesh = static_cast<int>(this->rcut/this->dr);
	if(this->mesh%2==0) this->mesh++;

	if(test==2)
	{
		cout << "\n ne=" << ne;
		cout << "\n lmax=" << lmax;
		cout << "\n mesh=" << mesh;
		cout << "\n total_nchi=" << total_nchi;
		cout << "\n dr=" << dr;
	}
	
	// init psi matrix, generate eigenvalue matrix
	// this eigenvalue will be used in many functions.
	this->eigenvalue.create( this->lmax+1, this->ne);
	this->jjnorm.create( this->lmax+1, this->ne);

	delete[] this->eigen1;
	this->eigen1 = new double[ne];

	delete[] this->r;
	delete[] this->rab;
	delete[] this->psi;
	delete[] this->psi_first;
	delete[] this->psi_second;
	delete[] this->below_max_os;
	delete[] this->os_number;

	this->r = new double[mesh];
	this->rab = new double[mesh];
	this->psi = new double[mesh];
	this->psi_first = new double[mesh];
	this->psi_second = new double[mesh];

	this->below_max_os = new bool[total_nchi];
	this->os_number = new int[total_nchi];
		
	for(int ir=0; ir<mesh; ir++)
	{
		r[ir] = dr * ir;
		rab[ir] = dr;
	}

	// (1) call to find the eigenvalues of Jlq.
	Calculate_C4::Find_Eigenvalues( 
			this->ne, 
			this->eigenvalue, 
			this->lmax, 
			this->rcut,
			tolerence);

	// (2) to calculate jl(qr)
	this->jle = new double**[this->lmax+1];
	double* g;
	if(SMOOTH) 
	{
		g = new double[mesh];
		double sigma2 = SIGMA*SIGMA;
		for(int ir=0; ir<mesh; ir++)
		{
			g[ir] = 1.0 - std::exp(-( (r[ir]-this->rcut)*(r[ir]-this->rcut)/2.0/sigma2) );
		}
	}
	
	for(int l=0; l<lmax+1; l++)
	{
		// get eigenvalues.
		for(int ie=0; ie<ne; ie++)
		{this->eigen1[ie] = this->eigenvalue(l, ie);}

		// allocate space
		this->jle[l] = new double*[ne];

		for(int ie=0; ie<this->ne; ie++)
		{
			this->jle[l][ie] = new double[this->mesh];
			Mathzone::Spherical_Bessel( this->mesh, this->r, this->eigen1[ie], l, jle[l][ie]);
			
			// (3) to calculate the jjnorm
			double* f = new double[mesh];
			for(int ir=0; ir<mesh; ir++)
			{
				f[ir] = std::pow(jle[l][ie][ir],2);
			}
			Mathzone::Simpson_Integral(mesh, f, rab, this->jjnorm(l,ie) );	
			jjnorm(l,ie)*=4.0*PI;
//			cout << "\n jj = " << jjnorm(l,ie);
			delete[] f;

			if(SMOOTH)
			for(int ir=0; ir<mesh; ir++)
			{
				jle[l][ie][ir] *= g[ir];
			}

		}
	}

	if(SMOOTH) delete[] g;

//	this->ofso.open("ORBITAL_OSCILATION.txt");
	this->ofsk.open("ORBITAL_KINETIC.txt");

	return;
}


double Psi_Second::get_ecut( double *c4, const int &L)
{
	timer::tick("Psi_Second","get_ecut");

	static double fpi = 4*PI;
	static double sfpi = sqrt(fpi);
	static double prefac_k = fpi / pow(2*PI,1.5);
	static double dk = 0.01;
	static int meshK = 800; // 8*8*2=128Ry
	static double* g;
	static double* rab;
	static double* rabk;
	static double* functionk;
	static bool init = false;
	static double norm;
	static double normk;
	static double ecut;

	if(!init)
	{
		if(meshK%2==0) ++meshK;
		
		g = new double[mesh]; // smooth function
		rab = new double[mesh];
		rabk = new double[meshK];
		functionk = new double[meshK];
		if(SMOOTH)
		{
			double sigma2 = SIGMA*SIGMA;
			for(int ir=0; ir<mesh; ir++)
			{
				g[ir] = 1.0 - std::exp(-( (r[ir]-this->rcut)*(r[ir]-this->rcut)/2.0/sigma2) );
			}
		}
		for(int ir=0; ir<mesh; ir++)
		{
			rab[ir] = dr;
		}
		for(int ik=0; ik<meshK; ik++)
		{
			rabk[ik] = dk;
		}
		init=true;
	}

	// get the eigenvalues.
	for(int ie=0; ie<ne; ie++)
	{this->eigen1[ie] = this->eigenvalue(L, ie);}

	ZEROS(psi,mesh);
	for(int ie=0; ie<ne; ie++)
	{
		for(int ir=0; ir<mesh; ir++)
		{
			psi[ir] += c4[ie] * jle[L][ie][ir];
		}
	}
	
	double *function = new double[mesh];
	for(int ir=0; ir<mesh; ir++)
	{
		function[ir] = pow(psi[ir]*r[ir],2);
	}

	Mathzone::Simpson_Integral(mesh, function, rab, norm );
//	cout << "\n norm = " << norm * fpi<< endl;
	assert(norm!=0.0);
//	exit(0);
	
	norm = sqrt(norm);
	for(int ir=0; ir<mesh; ir++)
	{
		psi[ir] /= norm;
	}

	for(int ie=0; ie<ne; ie++)
	{
		c4[ie] /= norm;
	}
	
	// NOTICE !!!! here psi ==> psi*r*r
	for(int ir=0; ir<mesh; ir++)
	{
		psi[ir] = psi[ir]/sfpi*r[ir]*r[ir];
	}

	// calculate psik.
	double *psik = new double[meshK];
	double *jj = new double[mesh];
	for(int ik=0; ik<meshK; ik++)
	{
		double q = ik*dk;
		Mathzone::Spherical_Bessel( mesh, r, q, L, jj);
		for(int ir=0; ir<mesh; ir++)
		{
			function[ir] = psi[ir]*jj[ir];
		} 
		Mathzone::Simpson_Integral(mesh, function, rab, psik[ik] );
	}	

	for(int ik=0; ik<meshK; ik++)
	{
		psik[ik] *= prefac_k;
	}

	//
	for(int ik=0; ik<meshK; ik++)
	{
		functionk[ik] = std::pow(psik[ik]*ik*dk,2);
	}
	Mathzone::Simpson_Integral(meshK, functionk, rabk, normk );
//	cout << "\n Psi(k) norm=" << normk * fpi;
	
	if(normk*fpi<0.999)
	{
		ecut = pow(dk*meshK,2)*2;	
	}
	else
	{
		int kmesh_used = meshK;
		while(normk*fpi> 0.999 && kmesh_used >= 0)
		{
			kmesh_used -= 2;
			Mathzone::Simpson_Integral(kmesh_used, functionk, rabk, normk );
		}
		ecut = pow(dk*kmesh_used,2)*2;
	}
	//cout << "Ecut(Ry)(norm>0.999)=" << ecut << endl;

	delete[] jj;
	delete[] psik;
	delete[] function;
	timer::tick("Psi_Second","get_ecut");
	return ecut;
}


void Psi_Second::norm_c4( double *c4, const int &L)
{
	ZEROS(this->psi,mesh);
	for(int ie=0; ie<ne; ie++)
	{
//		ZEROS(this->jle, mesh);
//		Mathzone::Spherical_Bessel( this->mesh, this->r, this->eigen1[ie], L, jle);
		
		for(int ir=0; ir<mesh; ir++)
		{this->psi[ir] += c4[ie] * jle[L][ie][ir];}
	}

	double* function = new double[mesh];
	for(int ir=0; ir<mesh; ir++)
	{
		function[ir] = pow(psi[ir]*r[ir],2);	
	}
	
	double norm;
	Mathzone::Simpson_Integral(mesh, function, rab, norm );	
	delete[] function;

	assert(norm>0.0);
	for(int ie=0; ie<ne; ie++)
	{
		c4[ie] /= std::sqrt(norm);
	}

	return;
}


