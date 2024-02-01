#include "ORB_atomic_lm.h"
#include "module_base/sph_bessel_recursive.h"
#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/constants.h"

#ifdef _OPENMP
#include <omp.h>
#endif

Numerical_Orbital_Lm::Numerical_Orbital_Lm()
{
	label = "";
	index_atom_type = 0;
	angular_momentum_l = 0;
	index_chi = 0;

	nr=1;
	nk=1;

	rcut=0.0;
	kcut=0.0;
	dk=0.0;

	nr_uniform = 1;
	dr_uniform = -1.0;
	zty = 0.0;	
}

Numerical_Orbital_Lm::~Numerical_Orbital_Lm()
{}

void Numerical_Orbital_Lm::set_orbital_info
(
 	const std::string &label_in,
	const int &index_atom_type_in,
	const int &angular_momentum_l_in,
	const int &index_chi_in,
	const int &nr_in,
	const double *rab_in,
	const double *r_radial_in,
	const Psi_Type &psi_type,				// Peize Lin add 2017-12-12
	const double *psi_in,
	const int &nk_in,
	const double &dk_in,
	// Peize Lin delete lat0 2016-02-03
	const double &dr_uniform_in,
	bool flag_plot,				// Peize Lin add flag_plot 2016-08-31
	bool flag_sbpool,			// Peize Lin add flag_sbpool 2017-10-02
	const bool &force_flag // mohan add 2021-05-07
)
{
	copy_parameter(
		label_in,
		index_atom_type_in,
		angular_momentum_l_in,
		index_chi_in,
		nr_in,
		rab_in,
		r_radial_in,
		nk_in,
		dk_in,
		dr_uniform_in);

	switch(psi_type)
	{
		case Psi_Type::Psi:
			for (int ir = 0; ir < nr; ir++)
			{
				this->psi[ir] = psi_in[ir];
				this->psir[ir] = psi[ir] * r_radial[ir]; //mohan 2010-04-19
			}
			break;
		case Psi_Type::Psif:
			for( int ik=0; ik!=nk; ++ik )
			{
				this->psif[ik] = psi_in[ik];
				this->psik[ik] = psif[ik] * k_radial[ik];
				this->psik2[ik] = psik[ik] * k_radial[ik];
			}
			break;
		case Psi_Type::Psik:
			psif.resize(0);
			for( int ik=0; ik!=nk; ++ik )
			{
				this->psik[ik] = psi_in[ik];
				this->psik2[ik] = psik[ik] * k_radial[ik];
			}
			break;
		case Psi_Type::Psik2:
			psif.resize(0);
			psik.resize(0);
			for( int ik=0; ik!=nk; ++ik )
				this->psik2[ik] = psi_in[ik];
			break;
		default:
			throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}

	switch(psi_type)
	{
		case Psi_Type::Psif:
		case Psi_Type::Psik:
		case Psi_Type::Psik2:
			if( flag_sbpool )
			{
			 	this->cal_rradial_sbpool();
			}
			else
			{
				throw std::domain_error("flag_sbpool false not finished in Numerical_Orbital_Lm::set_orbital_info_k. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
			}
			break;
		default:	break;
	}

	//liaochen modify on 2010/4/7
	//we do SBT on regular mesh
	//so we first generate psi_uniform first
	//we put uniform in ahead of cal_kradial
	
	/*
	bool uni = true;
	if (uni)
	{
		this->extra_uniform(dr_uniform, force_flag);
	}
	else
	{
		this->use_uniform(dr_uniform);
	}
	*/
	this->extra_uniform(dr_uniform, force_flag);

	switch(psi_type)
	{
		case Psi_Type::Psi:
			if( flag_sbpool )
			{
			 	this->cal_kradial_sbpool();
			}
			else
			{
				this->cal_kradial();
			}
			break;
		default:	break;
	}

//	this->norm_test();						// Peize Lin delete 2016-08-31
	if( flag_plot )
	{
		this->plot();			// Peize Lin add flag_plot 2016-08-31
	}
	return;
}

void Numerical_Orbital_Lm::copy_parameter(
 	const std::string &label_in,
	const int &index_atom_type_in,
	const int &angular_momentum_l_in,
	const int &index_chi_in,
	const int &nr_in,
	const double *rab_in,
	const double *r_radial_in,
	const int &nk_in,
	const double &dk_in,
	const double &dr_uniform_in)
{
    this->label = label_in;
    this->index_atom_type = index_atom_type_in;
    this->angular_momentum_l = angular_momentum_l_in;
    this->index_chi = index_chi_in;
    assert(nr_in>=2);
//  assert(nr_in<10000);									// Peize Lin delete 2017-12-03
    assert(nr%2!=0);
    this->nr = nr_in;
    assert(r_radial_in[nr-1]>0.0);
    // assert(r_radial_in[nr-1]<50);						// Peize Lin delete 2017-08-18
    this->rcut = r_radial_in[nr-1];
    assert(nk_in>1);
    //assert(nk_in<10000);				// Jiyy delete 2022-07-18
    this->nk = nk_in;
    assert(nk%2!=0);
    assert(dk_in>0);
    this->dk = dk_in;
	this->dr_uniform=dr_uniform_in;

	/***********************************************************
	be careful! LiaoChen modify on 2010/4/21
	************************************************************/
//	this->dk = ModuleBase::PI / rcut / 2.0;
//	this->nk = this->nr;

	r_radial.resize(nr);
	rab.resize(nr);
	psi.resize(nr);
	psir.resize(nr);
	for (int ir = 0; ir < nr; ir++)
	{
		this->r_radial[ir] = r_radial_in[ir];
		this->rab[ir] = rab_in[ir];
	}

	k_radial.resize(nk);
	psif.resize(nk);
	psik.resize(nk);
	psik2.resize(nk);
	for (int ik = 0; ik < nk; ik++)
	{
		this->k_radial[ik] = ik * this->dk;
	}
	this->kcut = (nk-1) * this->dk;
}

#include "module_base/mathzone_add1.h"
void Numerical_Orbital_Lm::extra_uniform(const double &dr_uniform_in, const bool &force_flag)
{
	ModuleBase::timer::tick("NOrbital_Lm", "extra_uniform");
	
	//---------------------------------------------
	// set the dr, fixed by liaochen.
	// calculate the number of radial mesh points.
	//---------------------------------------------
	assert(dr_uniform>0.0);
	this->dr_uniform = dr_uniform_in;
	this->nr_uniform = static_cast<int>(rcut/dr_uniform) + 10;
	
	this->psi_uniform.resize(nr_uniform,0);

	// do interpolation here to make grid more dense

#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (int ir = 0; ir < this->nr_uniform; ir++)
	{
		const double psi_uniform_tmp  = 
		ModuleBase::Mathzone_Add1::Uni_RadialF(ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->psi), this->nr, this->rab[0], ir * dr_uniform); 
		this->psi_uniform[ir] = psi_uniform_tmp;
//    	this->psi_uniform[ir] = ModuleBase::Mathzone::Polynomial_Interpolation(this->psi, this->nr, this->rab[0], ir * dr_uniform); 
    }
	
	//----------------------------------------------	 
	// calculate the dpsi_uniform
	//----------------------------------------------	 
	this->dpsi_uniform.resize(this->nr_uniform);
	this->ddpsi_uniform.resize(this->nr_uniform);

	double* y2 = new double[nr];

	//--------------------------------------------------------------------------
	// old code to calculate the derivate dpsi/dr, 
	// has problem that the derivatives of orbitals oscillate a lot
	// around r=0
	//--------------------------------------------------------------------------
	//ModuleBase::Mathzone_Add1::SplineD2 (r_radial, psi, nr, 100000.0, 100000.0, y2);
	//double yp1=(this->psi[1]-this->psi[0])/this->r_radial[1];
	//std::cout<<"psi0="<<"  "<<this->psi[0]<<"  "<<"psi1="<<"  "<<this->psi[1]<<"  "<<"r1="<<"  "<<this->r_radial[1]<<std::endl; 
	//std::cout<<"yp1="<<"  "<<yp1<<std::endl;
	//ModuleBase::Mathzone_Add1::SplineD2 (r_radial, psi, nr, yp1, 0.0, y2);
	

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// new code developed by pengfei.
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Peize Lin update 2016-08-31
	switch( this->angular_momentum_l ) // added by pengfei 13-8-8 different l has different  boundary conditions 
	{
		case 0: ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 0.0, 0.0, y2); break;
		case 1: ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 100000.0, 100000.0, y2); break;
		case 2: ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 0.0, 0.0, y2); break;
		case 3: ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 100000.0, 100000.0, y2); break;
		case 4: ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 0.0, 0.0, y2); break;
		default: 
			//GlobalV::ofs_warning << " The angular momentum larger than 4 (g orbitals) may be error about eggbox. " << std::endl;
			//GlobalV::ofs_warning << " Check file " << __FILE__ << " line " << __LINE__ <<std::endl;
			//std::cout << " The angular momentum larger than 4 (g orbitals) may be error about eggbox. " << std::endl;
			//std::cout << " Check file " << __FILE__ << " line " << __LINE__ <<std::endl;
			ModuleBase::Mathzone_Add1::SplineD2 (ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), nr, 0.0, 0.0, y2); break;
	}

	//ModuleBase::Mathzone_Add1::SplineD2 (r_radial, psi, nr, 0.0, 0.0, y2);
	//std::cout<<"angular_momentum_l="<<"  "<<this->angular_momentum_l<<std::endl;
	//for (int i=0; i<nr; i++)
	//{
	//     std::cout<<r_radial[i]<<"  "<<y2[i]<<std::endl;
	//}
	//Method 1
	//	ModuleBase::Mathzone_Add1::Uni_Deriv_Phi (psi_uniform, nr_uniform, dr_uniform, 1, dpsi_uniform);
	//	ModuleBase::Mathzone_Add1::Uni_Deriv_Phi (psi_uniform, nr_uniform, dr_uniform, 2, ddpsi_uniform);

	double* rad = new double[nr_uniform];
	for (int ir = 0; ir < nr_uniform; ir++)
	{
		rad[ir] = ir*dr_uniform;
	}

	//	ModuleBase::Mathzone_Add1::SplineD2 (rad, psi_uniform, nr_uniform, 0.0, 0.0, ddpsi_uniform);
	double* tmp = new double[nr_uniform];
	ModuleBase::Mathzone_Add1::Cubic_Spline_Interpolation(ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial), ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), y2, 
			nr, rad, nr_uniform, tmp, ModuleBase::GlobalFunc::VECTOR_TO_PTR(dpsi_uniform));

	// calculate zty
	// liaochen add 2010-08
	if( force_flag )	// Peize Lin add if 2017-10-26
	{
		ModuleBase::Mathzone_Add1::Uni_Deriv_Phi (
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->psi_uniform), 
			this->nr_uniform, 
			dr_uniform, 
			angular_momentum_l, 
			tmp);

		this->zty = tmp[0]/ModuleBase::Mathzone_Add1::factorial (angular_momentum_l);
	}

	delete [] y2;
	delete [] rad;
	delete [] tmp;
	ModuleBase::timer::tick("NOrbital_Lm", "extra_uniform");
}

/*
void Numerical_Orbital_Lm::use_uniform(const double &dr_uniform_in)
{
	assert(dr_uniform_in>0.0);
	this->dr_uniform = dr_uniform_in;
	// for save: +10, because in real space interpolation,
	// there may be "one grid point" more than the cutoff.
	this->nr_uniform = static_cast<int>(rcut/dr_uniform)+10;

	this->psi_uniform.resize(nr_uniform,0);

	std::string orbital_type; 
	// Peize Lin update 2016-08-31
	if( 0==this->angular_momentum_l )
	{
		orbital_type = 's';
	}
	else if( 1==this->angular_momentum_l )
	{
		orbital_type = 'p';
	}
	else if( 2==this->angular_momentum_l )
	{
		orbital_type = 'd';
	}
	else if( 3<=this->angular_momentum_l && this->angular_momentum_l<=6 )
	{
		orbital_type = 'f'+this->angular_momentum_l-3;
	}
	else if( 7<=this->angular_momentum_l && this->angular_momentum_l<=11 )
	{
		orbital_type = 'k'+this->angular_momentum_l-7;
	}
	else
	{
		orbital_type = "L" + ModuleBase::GlobalFunc::TO_STRING(this->angular_momentum_l);
	}
		
	std::cout << "===========================================================" << std::endl;
	for(int i=0; i<nr_uniform; i++)
	{
		this->psi_uniform[i] = 
			ModuleBase::Mathzone_Add1::Uni_RadialF(ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi), this->nr, this->rab[0], i*dr_uniform); 
	}
	
	this->dpsi_uniform.resize(nr_uniform);

	ModuleBase::Mathzone_Add1::Uni_Deriv_Phi (
		ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi_uniform), 
		nr_uniform, dr_uniform, 
		1, 
		ModuleBase::GlobalFunc::VECTOR_TO_PTR(dpsi_uniform));

#ifdef __NORMAL

#else	
	if(GlobalV::MY_RANK==0)
	{
		std::stringstream ss;
		ss << GlobalV::global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << ".ORBITAL_NOR_uniform.txt";

		std::ofstream ofs(ss.str().c_str());

		for(int i=0; i<nr_uniform; i++)
		{
			ofs << std::setw(15) << i*dr_uniform << std::setw(20) << psi_uniform[i] << std::endl;
		}
		ofs.close();
	}
#endif

	return;
}
*/

//liaochen modify on 2010/4/7
//use Sbt_new
void Numerical_Orbital_Lm::cal_kradial(void)
{
	assert( this->nr > 0);
	assert( this->nr_uniform > 0);
	double *jl = new double[nr];
	double *integrated_func = new double[nr];

	const double pref = sqrt( 2.0 / ModuleBase::PI );
	//Sbt method
	
	/*
	double* rad = new double[nr_uniform];
	for (int ir = 0; ir < nr_uniform; ir++) 
	{
		rad[ir] = dr_uniform * ir;
	}
	
	//liaochen add
	ModuleBase::Mathzone_Add1::Sbt_new (3, angular_momentum_l, 
							k_radial, dk, nk, 
							rad, dr_uniform, nr_uniform, 
							psi_uniform, 0, this->psik);
	
	for (int ik = 0; ik < nk; ik++) this->psik[ik] *= (pref*k_radial[ik]);
	delete [] rad;
	*/
	
	//integration directly
	for (int ik = 0; ik < nk; ik++)
	{
		ModuleBase::Sphbes::sphbesj(
				this->nr, 
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->r_radial), 
				this->k_radial[ik], 
				this->angular_momentum_l, 
				jl);

		for (int ir = 0; ir < nr; ir++)
		{
			integrated_func[ir] = this->psir[ir] * this->r_radial[ir] * jl[ir];
		}

		ModuleBase::Integral::Simpson_Integral(
				this->nr,
				integrated_func,
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->rab),
				this->psif[ik]);
		this->psif[ik] *= pref;
		this->psik[ik] = this->psif[ik] * k_radial[ik];
		this->psik2[ik] = this->psik[ik] * k_radial[ik];
	}
		
	delete[] integrated_func;
	delete[] jl;
}

/*
// Peize Lin add 2017-10-02
void Numerical_Orbital_Lm::cal_kradial_sbpool(void)
{
	assert( this->nr > 0);
	assert( this->nr_uniform > 0);

	// dr must be all the same for Sph_Bessel_Recursive_Pool and Simpson_Integral
	const double dr = this->rab[0];
	for( size_t ir=1; ir<this->nr; ++ir )
		assert( dr == this->rab[ir] );

	ModuleBase::Sph_Bessel_Recursive::D2* pSB = nullptr;
	for( auto & sb : Sph_Bessel_Recursive_Pool::D2::sb_pool )
		if( this->dk * dr == sb.get_dx() )
		{
			pSB = &sb;
			break;
		}
	if(!pSB)
	{
		Sph_Bessel_Recursive_Pool::D2::sb_pool.push_back({});
		pSB = &Sph_Bessel_Recursive_Pool::D2::sb_pool.back();
	}
	pSB->set_dx( this->dk * dr );
	pSB->cal_jlx( this->angular_momentum_l, this->nk, this->nr );
	const std::vector<std::vector<double>> &jl = pSB->get_jlx()[this->angular_momentum_l];

	std::vector<double> integrated_func( this->nr );
	const double pref = sqrt( 2.0 / ModuleBase::PI );

	std::vector<double> psir2(nr);
	for( size_t ir=0; ir!=nr; ++ir )
		psir2[ir] = this->psir[ir] * this->r_radial[ir];

	for (int ik = 0; ik < nk; ik++)
	{
		const std::vector<double> &jlk = jl[ik];
		for (int ir = 0; ir < nr; ir++)
			integrated_func[ir] = psir2[ir] * jlk[ir];
		ModuleBase::Integral::Simpson_Integral(
				this->nr,
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func),
				dr,
				this->psik[ik]);
		this->psik[ik] *= ( pref * k_radial[ik]);
	}
}
*/

// Peize Lin add 2017-10-27
void Numerical_Orbital_Lm::cal_kradial_sbpool(void)
{
	assert( this->nr > 0);
	assert( this->nr_uniform > 0);

	// dr must be all the same for Sph_Bessel_Recursive_Pool
	const double dr = this->rab[0];
	
	for( int ir=1; ir<this->nr; ++ir )
	{
		assert( dr == this->rab[ir] );
	}

	ModuleBase::Sph_Bessel_Recursive::D2* pSB = nullptr;
	for( auto & sb : ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool )
	{
		if( this->dk * dr == sb.get_dx() )
		{
			pSB = &sb;
			break;
		}
	}

	if(!pSB)
	{
		ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.push_back({});
		pSB = &ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.back();
	}
	pSB->set_dx( this->dk * dr );
	pSB->cal_jlx( this->angular_momentum_l, this->nk, this->nr );
	const std::vector<std::vector<double>> &jl = pSB->get_jlx()[this->angular_momentum_l];

	const double pref = sqrt( 2.0 / ModuleBase::PI );

	std::vector<double> r_tmp(nr);
	for( int ir=0; ir!=nr; ++ir )
	{
		r_tmp[ir] = this->psir[ir] * this->r_radial[ir] * this->rab[ir];
	}

	constexpr double one_three=1.0/3.0, two_three=2.0/3.0, four_three=4.0/3.0;
	r_tmp[0]*=one_three;	
	r_tmp[nr-1]*=one_three;

	for( int ir=1; ir!=nr-1; ++ir )
	{
		r_tmp[ir] *= (ir&1) ? four_three : two_three;
	}
#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (int ik = 0; ik < nk; ik++)
	{
#ifdef __NORMAL
		double psi_f_tmp = 0.0; 
		for(int ir=0; ir<nr; ++ir)
		{
			psi_f_tmp += r_tmp[ir]*jl[ik][ir];
		}
		psi_f_tmp *= pref;
#else
		const double psi_f_tmp = 
		pref * BlasConnector::dot( this->nr, ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_tmp), 1, ModuleBase::GlobalFunc::VECTOR_TO_PTR(jl[ik]), 1 ) ;
#endif
		this->psif[ik] = psi_f_tmp;
		this->psik[ik] = psi_f_tmp * k_radial[ik];
		this->psik2[ik] = this->psik[ik] * k_radial[ik];
	}
	return;
}

// Peize Lin add 2017-12-11
void Numerical_Orbital_Lm::cal_rradial_sbpool(void)
{
	// dr must be all the same for Sph_Bessel_Recursive_Pool
	const double dr = this->rab[0];

	for( int ir=1; ir<this->nr; ++ir )
	{
		assert( dr == this->rab[ir] );
	}

	ModuleBase::Sph_Bessel_Recursive::D2* pSB = nullptr;
	for( auto & sb : ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool )
	{
		if( dr * dk == sb.get_dx() )
		{
			pSB = &sb;
			break;
		}
	}

	if(!pSB)
	{
		ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.push_back({});
		pSB = &ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.back();
	}

	pSB->set_dx( dr * dk );
	pSB->cal_jlx( this->angular_momentum_l, this->nr, this->nk );

	const std::vector<std::vector<double>> &jl = pSB->get_jlx()[this->angular_momentum_l];

	const double pref = sqrt(2.0/ModuleBase::PI);

	std::vector<double> k_tmp(nk);

	for( int ik=0; ik!=nk; ++ik )
	{
		k_tmp[ik] = this->psik2[ik] * dk;
	}

	constexpr double one_three=1.0/3.0, two_three=2.0/3.0, four_three=4.0/3.0;

	k_tmp[0]*=one_three;	
	k_tmp[nk-1]*=one_three;

	for( int ik=1; ik!=nk-1; ++ik )
	{
		k_tmp[ik] *= (ik&1) ? four_three : two_three;
	}

	for( int ir = 0; ir!=nr; ++ir )
	{
#ifdef __NORMAL
		// mohan add 2021-05-08, test needed
		double kj_dot = 0.0;
		for( int ik=0; ik<nk; ++ik)
		{
			kj_dot += k_tmp[ik]*jl[ir][ik];
		}
		this->psi[ir] = pref * kj_dot;
#else
		this->psi[ir] = pref * BlasConnector::dot( this->nk, ModuleBase::GlobalFunc::VECTOR_TO_PTR(k_tmp), 1, ModuleBase::GlobalFunc::VECTOR_TO_PTR(jl[ir]), 1 );
#endif
		this->psir[ir] = this->psi[ir] * r_radial[ir];
	}
}

//===============================================
//FOUND LOCAL VARIABLE
//asum : integral of psi*psi in whole space
//===============================================
/*
void Numerical_Orbital_Lm::norm_test(void)const
{
//	ModuleBase::TITLE(ofs_onscaling, "Numerical_Orbital_Lm", "norm_test");
	//double asum_r = 0.0;
	//double asum_k = 0.0;

	// note here psir = psi * r
	double *f = new double[nr];
	for(int ir=0; ir<nr; ir++)
	{
		f[ir] = this->psir[ir] * this->psir[ir];
	}

	double sumr = 0.0;
	//double sumk = 0.0;

	ModuleBase::Integral::Simpson_Integral(this->nr, f, ModuleBase::GlobalFunc::VECTOR_TO_PTR(this->rab), sumr);

	delete[] f;
	f = new double[nk];
	for(int ik=0; ik<nk; ik++)
	{
		f[ik] = this->psik[ik] * this->psik[ik];
	}

//	ModuleBase::Integral::Simpson_Integral(this->nk, f, this->k_radial, sumk);
	
	//means nothing.
	//GlobalV::ofs_running << std::setw(12) << sumk << std::endl;

	delete[] f;
	return;
}
*/

void Numerical_Orbital_Lm::plot(void)const
{
	ModuleBase::TITLE("Numerical_Orbital_Lm","plot");
	
	std::string orbital_type;
	// Peize Lin update 2016-08-31
	if( 0==this->angular_momentum_l )
	{
		orbital_type = 's';
	}
	else if( 1==this->angular_momentum_l )
	{
		orbital_type = 'p';
	}
	else if( 2==this->angular_momentum_l )
	{
		orbital_type = 'd';
	}
	else if( 3<=this->angular_momentum_l && this->angular_momentum_l<=6 )
	{
		orbital_type = 'f' + this->angular_momentum_l - 3;
	}
	else if( 7<=this->angular_momentum_l && this->angular_momentum_l<=11 )
	{
		orbital_type = 'k' + this->angular_momentum_l - 7;
	}
	else
	{
		orbital_type = "L" + ModuleBase::GlobalFunc::TO_STRING(this->angular_momentum_l);	
	}

	if(GlobalV::MY_RANK==0)
	{
		std::stringstream ssr, ssk, ssru ,ssdru; // 2013-08-10 pengfei
		ssr << GlobalV::global_out_dir << this->label << "/"
			<< this->label << "-"<< orbital_type << index_chi+1 << "-orbital-r.dat";

		ssk << GlobalV::global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-k.dat";

		ssru << GlobalV::global_out_dir << this->label << "/"
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-ru.dat";

		ssdru << GlobalV::global_out_dir << this->label << "/"  // 2013-08-10 pengfei
			<< this->label << "-" << orbital_type << index_chi+1 << "-orbital-dru.dat";

		std::ofstream ofsr(ssr.str().c_str());
		std::ofstream ofsk(ssk.str().c_str());
		std::ofstream ofsru(ssru.str().c_str());
		std::ofstream ofsdru(ssdru.str().c_str()); // 2013-08-10 pengfei

		if (!ofsk || !ofsr || !ofsru || !ofsdru) // 2013-08-10 pengfei
		{
			ModuleBase::WARNING("Numerical_Orbital_Lm : plot", "Can't open files !");
		}

		for (int i = 0; i < this->nr; i++)
		{
			ofsr << this->r_radial[i] << " " << psi[i] << std::endl;
		}

		for (int i = 0; i < this->nk; i++)
		{
			ofsk << this->k_radial[i] << " " << psik[i] << std::endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsru << this->dr_uniform * i << " " << psi_uniform[i] << std::endl;
		}

		for (int i = 0; i < this->nr_uniform; i++)
		{
			ofsdru << this->dr_uniform * i << " " << dpsi_uniform[i] << std::endl;// output dphi/dr 2013-08-10  pengfei
		}
		ofsr.close();
		ofsk.close();
		ofsru.close();
		ofsdru.close(); // 13-08-10 pengfei
	}

	return;
}
