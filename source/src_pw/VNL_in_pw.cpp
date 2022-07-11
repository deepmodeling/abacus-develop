#include "global.h"
#include "VNL_in_pw.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "wavefunc.h"
#include "../module_orbital/ORB_gen_tables.h"
#include "../module_base/math_integral.h"
#include "../module_base/math_sphbes.h"
#include "../module_base/math_polyint.h"
#include "../module_base/math_ylmreal.h"
#include "soc.h"
#include "../module_base/timer.h"

pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
#ifdef __CUDA
	cudaFree(this->d_deeq);
#endif
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq 
// allocate vkb, GlobalV::NQX, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const int ntype, const bool allocate_vkb)
{
	ModuleBase::TITLE("pseudopot_cell_vnl", "init");
	ModuleBase::timer::tick("ppcell_vnl", "init");

	GlobalV::ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << std::endl;

	int it = 0;
//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : lmaxkb(max angular momentum,(see pseudo_h))
//----------------------------------------------------------
	this->lmaxkb = - 1;
	for (it = 0;it < ntype; it++)
	{
		GlobalV::ofs_running << " " << GlobalC::ucell.atoms[it].label << " non-local projectors:" << std::endl;
		for (int ibeta = 0; ibeta < GlobalC::ucell.atoms[it].nbeta; ibeta++) 
		{
			GlobalV::ofs_running << " projector " << ibeta+1 << " L=" << GlobalC::ucell.atoms[it].lll[ibeta] <<  std::endl;
			this->lmaxkb = std::max( this->lmaxkb, GlobalC::ucell.atoms[it].lll[ibeta]);
		}
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nhm(max number of different beta functions per atom)
//----------------------------------------------------------
	this->nhm = 0;
	for (it=0;it<ntype;it++)
	{	
		this->nhm = std::max(nhm, GlobalC::ucell.atoms[it].nh);
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nkb(total number of beta functions, with struct.fact.)
//----------------------------------------------------------
	this->nkb = 0;
	for (it=0; it<ntype; it++)
	{
		this->nkb += GlobalC::ucell.atoms[it].nh * GlobalC::ucell.atoms[it].na;
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"TOTAL NUMBER OF NONLOCAL PROJECTORS",nkb);

	if( this->nhm > 0 )
	{
		this->indv.create(ntype, this->nhm);
		this->nhtol.create(ntype, this->nhm);
		this->nhtolm.create(ntype, this->nhm);
		this->nhtoj.create(ntype, this->nhm);
		this->deeq.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
#ifdef __CUDA
		cudaMalloc((void**)&d_deeq, GlobalV::NSPIN*GlobalC::ucell.nat*this->nhm*this->nhm*sizeof(double));
#endif		
		this->deeq_nc.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
		this->dvan.create(ntype, this->nhm, this->nhm);
		this->dvan_so.create(GlobalV::NSPIN, ntype, this->nhm, this->nhm);
		this->becsum.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm * (this->nhm + 1) / 2);
	}
	else
	{
		GlobalV::ofs_running << "\n nhm = 0, not allocate some matrix.";
	}

	// nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
	// dq+4)*cell_factor;
	this->lmaxq = 2 * this->lmaxkb + 1;

	if (nkb > 0 && allocate_vkb )
	{
		vkb.create(nkb, GlobalC::wf.npwx);
	}

	//this->nqx = 10000;		// calculted in allocate_nlpot.f90
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ); //LiuXh modify 20180515
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) + 1000; //LiuXh add 20180515
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) * 10; //LiuXh add 20180515
	GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) * cell_factor; //LiuXh add 20180619
	// nqx = (sqrt(ecutwfc)/dq+4)*cell_factor;

	
	// mohan update 2021-02-22
	const int nbrx = 10;
	const int nbrx_nc = 20;
	//  max number of beta functions
	if(GlobalV::NSPIN!=4) 
	{
		this->tab.create(ntype, nbrx, GlobalV::NQX);
	}
	else 
	{
		this->tab.create(ntype, nbrx_nc, GlobalV::NQX);
	}

	
	// mohan update 2021-02-22
	int nchix = 10;
	int nchix_nc = 20;
	// nchix : max number of atomic wavefunctions per atom
	if(GlobalV::NSPIN!=4) 
	{
		this->tab_at.create(ntype, nchix, GlobalV::NQX);
	}
	else 
	{
		this->tab_at.create(ntype, nchix_nc, GlobalV::NQX);
	}

	ModuleBase::timer::tick("ppcell_vnl","init");
	return;
}



//----------------------------------------------------------
// Calculates beta functions (Kleinman-Bylander projectors),
// with structure factor, for all atoms, in reciprocal space
//----------------------------------------------------------
void pseudopot_cell_vnl::getvnl(const int &ik, ModuleBase::ComplexMatrix& vkb_in)const
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl");
	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	if(lmaxkb < 0) 
	{
		return;
	}

	const int npw = GlobalC::kv.ngk[ik];

	ModuleBase::matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix ylm(x1, npw);
	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (int ig = 0;ig < npw;ig++) 
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}

	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].nbeta;
		const int nh = GlobalC::ucell.atoms[it].nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (int nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (int ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;

				vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
						this->tab, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );
			}

			// add spherical harmonic part
			for (int ih = 0;ih < nh;ih++)
			{
				if (nb == this->indv(it, ih))
				{
					const int lm = static_cast<int>( nhtolm(it, ih) );
					for (int ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++) 
		{
			std::complex<double> *sk = GlobalC::wf.get_sk(ik, it, ia,GlobalC::wfcpw);
			for (int ih = 0;ih < nh;ih++)
			{
				std::complex<double> pref = pow( ModuleBase::NEG_IMAG_UNIT, nhtol(it, ih));	//?
				std::complex<double>* pvkb = &vkb_in(jkb, 0);
				for (int ig = 0;ig < npw;ig++)
				{
					pvkb[ig] = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;

	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	return;
} // end subroutine getvnl



void pseudopot_cell_vnl::init_vnl(UnitCell_pseudo &cell)
{
	ModuleBase::TITLE("pseudopot_cell_vnl","init_vnl");
	ModuleBase::timer::tick("ppcell_vnl","init_vnl");

	//from init_us_1
	//   a) For each non vanderbilt pseudopotential it computes the D and
	//      the betar in the same form of the Vanderbilt pseudopotential.
	//   b) It computes the indices indv which establish the correspondence
	//      nh <-> beta in the atom
	//   c) It computes the indices nhtol which establish the correspondence
	//      nh <-> angular momentum of the beta function
	//   d) It computes the indices nhtolm which establish the correspondence
	//      nh <-> combined (l,m) index for the beta function.

	// For each pseudopotential we initialize the indices nhtol, nhtolm,
	// nhtoj, indv, and if the pseudopotential is of KB type we initialize
	// the atomic D terms

	this->dvan.zero_out();
	this->dvan_so.zero_out();//added by zhengdy-soc

	for(int it=0;it<cell.ntype;it++)
	{
		int BetaIndex=0;
		const int Nprojectors = cell.atoms[it].nh;
		for (int ib=0; ib<cell.atoms[it].nbeta; ib++)
		{
			const int l = cell.atoms[it].lll [ib];
			const double j = cell.atoms[it].jjj [ib];
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
				this->nhtoj(it,BetaIndex) = j;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}

		//    From now on the only difference between KB and US pseudopotentials
		//    is in the presence of the q and Q functions.
		//    Here we initialize the D of the solid
		if(cell.atoms[it].has_so )
		{
			Soc soc;
			soc.rot_ylm(this->lmaxkb);
			soc.fcoef.create(cell.ntype, this->nhm, this->nhm);
			for(int ip=0; ip<Nprojectors; ip++)
			{
				const int l1 = this->nhtol (it, ip);
				const double j1 = this->nhtoj (it, ip);
				const int m1 = this->nhtolm (it, ip) - l1* l1;
				//const int v1 = static_cast<int>( indv(it, ip ) );
				for(int ip2=0;ip2<Nprojectors; ip2++)
				{
					const int l2 = this->nhtol (it, ip2);
					const double j2 = this->nhtoj (it, ip2);
					const int m2 = this->nhtolm (it, ip2) - l2* l2;
					//const int v2 = static_cast<int>( indv(it, ip2 ) );
					if(l1 == l2 && fabs(j1-j2)<1e-7)
					{
						for(int is1=0;is1<2;is1++)
						{
							for(int is2=0;is2<2;is2++)
							{
								soc.set_fcoef(l1, l2,
										is1, is2,
										m1, m2,
										j1, j2,
										it, ip, ip2);
							}
						}
					}
				}
			}
//
//   and calculate the bare coefficients
//
			for(int ip = 0;ip<Nprojectors; ++ip)
			{
				const int ir = static_cast<int>( indv(it, ip ) );
				for(int ip2=0; ip2<Nprojectors; ++ip2)
				{
					const int is = static_cast<int>( indv(it, ip2) );
					int ijs =0;
					for(int is1=0;is1<2;++is1)
					{
						for(int is2=0;is2<2;++is2)
						{
							this->dvan_so(ijs,it,ip,ip2) = cell.atoms[it].dion(ir, is) * soc.fcoef(it,is1,is2,ip,ip2);
							++ijs;
							if(ir != is) soc.fcoef(it,is1,is2,ip,ip2) = std::complex<double>(0.0,0.0);
						}
					}
				}
			}
		}
		else
		for (int ip=0; ip<Nprojectors; ip++)
		{
			for (int ip2=0; ip2<Nprojectors; ip2++)
			{
				if ( this->nhtol (it, ip) == nhtol (it, ip2) &&
				     this->nhtolm(it, ip) == nhtolm(it, ip2) )
				{
					const int ir = static_cast<int>( indv(it, ip ) );
					const int is = static_cast<int>( indv(it, ip2) );
					if(GlobalV::LSPINORB)
					{
						this->dvan_so(0,it,ip,ip2) = cell.atoms[it].dion(ir, is);
						this->dvan_so(3,it,ip,ip2) = cell.atoms[it].dion(ir, is);
					}
					else
					{
						this->dvan(it, ip, ip2) = cell.atoms[it].dion(ir, is);
					}
				}
			} 
		} 
	} 

	// h) It fills the interpolation table for the beta functions
	/**********************************************************
	// He Lixin: this block is used for non-local potential
	// fill the interpolation table tab
	************************************************************/

	const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
	this->tab.zero_out();
	GlobalV::ofs_running<<"\n Init Non-Local PseudoPotential table : ";
	for (int it = 0;it < cell.ntype;it++)  
	{
		const int nbeta = cell.atoms[it].nbeta;
		int kkbeta = cell.atoms[it].kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			const int l = cell.atoms[it].lll[ib];
			for (int iq=0; iq<GlobalV::NQX; iq++)  
			{
				const double q = iq * GlobalV::DQ;
				ModuleBase::Sphbes::Spherical_Bessel(kkbeta, cell.atoms[it].r, q, l, jl);

				for (int ir = 0;ir < kkbeta;ir++)
				{
					aux[ir] = cell.atoms[it].betar(ib, ir) *
					          jl[ir] * cell.atoms[it].r[ir];
				} 
				double vqint;
				ModuleBase::Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].rab, vqint);
				this->tab(it, ib, iq) = vqint * pref;
			} 
		} 
		delete[] aux;
		delete[] jl;
	}
	ModuleBase::timer::tick("ppcell_vnl","init_vnl");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done." << std::endl;
	return;
}

#ifdef __LCAO
std::complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M)   // pengfei Li  2018-3-23
{
	std::complex<double> cf;
	if(alpha == 0)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,1,L,M);
	}
	else if(alpha == 1)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,2,L,M);
	}
	else if(alpha == 2)
	{
		cf = sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,0,L,M);
	}
	else
	{
		ModuleBase::WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
	}
	
	return cf;
}

double pseudopot_cell_vnl::CG(int l1, int m1, int l2, int m2, int L, int M)      // pengfei Li 2018-3-23
{
	int dim = L*L+M;
	int dim1 = l1*l1+m1;
	int dim2 = l2*l2+m2;
	
	//double A = MGT.Gaunt_Coefficients(dim1, dim2, dim);
	
	return MGT.Gaunt_Coefficients(dim1, dim2, dim);
}

void pseudopot_cell_vnl::getvnl_alpha(const int &ik)           // pengfei Li  2018-3-23
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl_alpha");
	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");

	if(lmaxkb < 0) 
	{
		return;
	}
	
	const int npw = GlobalC::kv.ngk[ik];
	int ig, ia, nb, ih, lu, mu;

	double *vq = new double[npw];
	const int x1= (lmaxkb + 2)*(lmaxkb + 2);

	ModuleBase::matrix ylm(x1, npw);
	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++) 
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}

	vkb1_alpha = new std::complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		vkb1_alpha[i] = new std::complex<double>*[nhm];
		for(int j=0; j<nhm; j++)
		{
			vkb1_alpha[i][j] = new std::complex<double>[npw];
		}
	}	
	
	vkb_alpha = new std::complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		vkb_alpha[i] = new std::complex<double>*[nkb];
		for(int j=0; j<nkb; j++)
		{
			vkb_alpha[i][j] = new std::complex<double>[GlobalC::wf.npwx];
		}
	}
	
	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

	MGT.init_Gaunt_CH( lmaxkb + 2 );
	MGT.init_Gaunt( lmaxkb + 2 );
	
	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].nbeta;
		const int nh = GlobalC::ucell.atoms[it].nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for(int i=0; i<3; i++)
			for(int j=0; j<nhm; j++)
			{
				ModuleBase::GlobalFunc::ZEROS(vkb1_alpha[i][j], npw);
			}
			
		for (ih = 0;ih < nh; ih++)
		{
			lu = static_cast<int>( nhtol(it, ih));
			mu = static_cast<int>( nhtolm(it, ih)) - lu * lu;
			nb = static_cast<int>( indv(it, ih));
			
			for (int L= abs(lu - 1); L<= (lu + 1); L++)
			{
				for (ig = 0;ig < npw;ig++)
				{
					const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
					vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
							this->tab_alpha, it, nb, L, GlobalV::NQX, GlobalV::DQ, gnorm);
					
					for (int M=0; M<2*L+1; M++)
					{
						int lm = L*L + M;
						for (int alpha=0; alpha<3; alpha++)
						{
							std::complex<double> c = Cal_C(alpha,lu, mu, L, M);
							/*if(alpha == 0)
							{
								std::cout<<"lu= "<<lu<<"  mu= "<<mu<<"  L= "<<L<<"  M= "<<M<<" alpha = "<<alpha<<"  "<<c<<std::endl;
							}*/
							vkb1_alpha[alpha][ih][ig] += c * vq[ig] * ylm(lm, ig) * pow( ModuleBase::NEG_IMAG_UNIT, L);
						}	
					}
				}
			}
		} // end nbeta

		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++) 
		{
			std::complex<double> *sk = GlobalC::wf.get_sk(ik, it, ia,GlobalC::wfcpw);
			for (ih = 0;ih < nh;ih++)
			{
				for (ig = 0;ig < npw;ig++)
				{
					for(int alpha=0; alpha<3; alpha++)
					{
						vkb_alpha[alpha][jkb][ig] = vkb1_alpha[alpha][ih][ig] * sk [ig];
					}					
				}
				++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;
	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");
	return;
} 
#endif

void pseudopot_cell_vnl::init_vnl_alpha(void)          // pengfei Li 2018-3-23
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","init_vnl_alpha");
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");

	for(int it=0;it<GlobalC::ucell.ntype;it++)
	{
		int BetaIndex=0;
		//const int Nprojectors = GlobalC::ucell.atoms[it].nh;
		for (int ib=0; ib<GlobalC::ucell.atoms[it].nbeta; ib++)
		{
			const int l = GlobalC::ucell.atoms[it].lll [ib];
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}
	} 


	// max number of beta functions
	const int nbrx = 10;

	const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);
	this->tab_alpha.create(GlobalC::ucell.ntype, nbrx, lmaxkb+2, GlobalV::NQX);
	this->tab_alpha.zero_out();
	GlobalV::ofs_running<<"\n Init Non-Local PseudoPotential table( including L index) : ";
	for (int it = 0;it < GlobalC::ucell.ntype;it++)  
	{
		const int nbeta = GlobalC::ucell.atoms[it].nbeta;
		int kkbeta = GlobalC::ucell.atoms[it].kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			for (int L = 0; L <= lmaxkb+1; L++)
			{
				for (int iq = 0; iq < GlobalV::NQX; iq++)
				{
					const double q = iq * GlobalV::DQ;
					ModuleBase::Sphbes::Spherical_Bessel(kkbeta, GlobalC::ucell.atoms[it].r, q, L, jl);
					
					for (int ir = 0;ir < kkbeta;ir++)
					{
						aux[ir] = GlobalC::ucell.atoms[it].betar(ib, ir) * jl[ir] * 
								  GlobalC::ucell.atoms[it].r[ir] * GlobalC::ucell.atoms[it].r[ir];
					}
					double vqint;
					ModuleBase::Integral::Simpson_Integral(kkbeta, aux, GlobalC::ucell.atoms[it].rab, vqint);
					this->tab_alpha(it, ib, L, iq) = vqint * pref;
				}
			}
		} 
		delete[] aux;
		delete[] jl;
	}
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << std::endl;
	return;
}



void pseudopot_cell_vnl::print_vnl(std::ofstream &ofs)
{
	output::printr3_d(ofs, " tab : ", tab);
}



int pseudopot_cell_vnl::calculate_nqx(const double &ecutwfc,const double &dq)
{
	int points_of_table = static_cast<int>( sqrt(ecutwfc)/dq + 4 ) ;//* cell_factor;
//	std::cout<<"\n nqx = "<< points_of_table;
//----------------------------------------------------------
// EXPLAIN : Plus 1 because the formula is transfered from 
// fortran code.
//----------------------------------------------------------
	return points_of_table + 1; 
}
