#include "global.h"
#include "VNL_in_pw.h"
#include "tools.h"
#include "wavefunc.h"
#include "../src_lcao/ORB_gen_tables.h"
#include "../src_global/math_integral.h"

pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq 
// allocate vkb, NQX, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const int ntype, const bool allocate_vkb)
{
	TITLE("pseudopot_cell_vnl", "init");
	timer::tick("ppcell_vnl","init",'C');

	ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << endl;

	int it = 0;
//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : lmaxkb(max angular momentum,(see pseudo_h))
//----------------------------------------------------------
	this->lmaxkb = - 1;
	for (it = 0;it < ntype; it++)
	{
		ofs_running << " " << ucell.atoms[it].label << " non-local projectors:" << endl;
		for (int ibeta = 0; ibeta < ucell.atoms[it].nbeta; ibeta++) 
		{
			ofs_running << " projector " << ibeta+1 << " L=" << ucell.atoms[it].lll[ibeta] <<  endl;
			this->lmaxkb = max( this->lmaxkb, ucell.atoms[it].lll[ibeta]);
		}
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nhm(max number of different beta functions per atom)
//----------------------------------------------------------
	this->nhm = 0;
	for (it=0;it<ntype;it++)
	{	
		this->nhm = std::max(nhm, ucell.atoms[it].nh);
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nkb(total number of beta functions, with struct.fact.)
//----------------------------------------------------------
	this->nkb = 0;
	for (it=0; it<ntype; it++)
	{
		this->nkb += ucell.atoms[it].nh * ucell.atoms[it].na;
	}

	OUT(ofs_running,"TOTAL NUMBER OF NONLOCAL PROJECTORS",nkb);

	if( this->nhm > 0 )
	{
		this->indv.create(ntype, this->nhm);
		this->nhtol.create(ntype, this->nhm);
		this->nhtolm.create(ntype, this->nhm);
		this->nhtoj.create(ntype, this->nhm);
		this->deeq.create(NSPIN, ucell.nat, this->nhm, this->nhm);
		this->deeq_nc.create(NSPIN, ucell.nat, this->nhm, this->nhm);
		this->dvan.create(ntype, this->nhm, this->nhm);
		this->dvan_so.create(NSPIN, ntype, this->nhm, this->nhm);
		this->becsum.create(NSPIN, ucell.nat, this->nhm * (this->nhm + 1) / 2);
	}
	else
	{
		ofs_running << "\n nhm = 0, not allocate some matrix.";
	}

	// nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
	// dq+4)*cell_factor;
	this->lmaxq = 2 * this->lmaxkb + 1;

	if (nkb > 0 && allocate_vkb )
	{
		vkb.create(nkb, wf.npwx);
	}

	//this->nqx = 10000;		// calculted in allocate_nlpot.f90
	//NQX = this->calculate_nqx(pw.ecutwfc,DQ); //LiuXh modify 20180515
	//NQX = this->calculate_nqx(pw.ecutwfc,DQ) + 1000; //LiuXh add 20180515
	//NQX = this->calculate_nqx(pw.ecutwfc,DQ) * 10; //LiuXh add 20180515
	NQX = this->calculate_nqx(pw.ecutwfc,DQ) * cell_factor; //LiuXh add 20180619
	// nqx = (sqrt(ecutwfc)/dq+4)*cell_factor;

	
	// mohan update 2021-02-22
	const int nbrx = 10;
	const int nbrx_nc = 20;
	//  max number of beta functions
	if(NSPIN!=4) 
	{
		this->tab.create(ntype, nbrx, NQX);
	}
	else 
	{
		this->tab.create(ntype, nbrx_nc, NQX);
	}

	
	// mohan update 2021-02-22
	int nchix = 10;
	int nchix_nc = 20;
	// nchix : max number of atomic wavefunctions per atom
	if(NSPIN!=4) 
	{
		this->tab_at.create(ntype, nchix, NQX);
	}
	else 
	{
		this->tab_at.create(ntype, nchix_nc, NQX);
	}

	timer::tick("ppcell_vnl","init",'C');
	return;
}



//----------------------------------------------------------
// Calculates beta functions (Kleinman-Bylander projectors),
// with structure factor, for all atoms, in reciprocal space
//----------------------------------------------------------
void pseudopot_cell_vnl::getvnl(const int &ik)
{
	if(test_pp) TITLE("pseudopot_cell_vnl","getvnl");
	timer::tick("pp_cell_vnl","getvnl");

	if(lmaxkb < 0) 
	{
		return;
	}

	const int npw = kv.ngk[ik];
	int ig, ia, nb, ih;
	matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	matrix ylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++) 
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}

	Mathzone::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		// calculate beta in G-space using an interpolation table
		const int nbeta = ucell.atoms[it].nbeta;
		const int nh = ucell.atoms[it].nh;

		if(test_pp>1) OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(test_pp>1) OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * ucell.tpiba;

				vq [ig] = Mathzone::Polynomial_Interpolation(
						this->tab, it, nb, NQX, DQ, gnorm );
			}

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == this->indv(it, ih))
				{
					const int lm = static_cast<int>( nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<ucell.atoms[it].na; ia++) 
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				complex<double> pref = pow( NEG_IMAG_UNIT, nhtol(it, ih));	//?
				for (ig = 0;ig < npw;ig++)
				{
					this->vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;

	timer::tick("pp_cell_vnl","getvnl");

	return;
} // end subroutine getvnl



void pseudopot_cell_vnl::init_vnl(UnitCell_pseudo &cell)
{
	TITLE("pseudopot_cell_vnl","init_vnl");
	timer::tick("ppcell_vnl","init_vnl");

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
	soc.rot_ylm(this->lmaxkb);
	soc.fcoef.create(cell.ntype, this->nhm, this->nhm);

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
							if(ir != is) soc.fcoef(it,is1,is2,ip,ip2) = complex<double>(0.0,0.0);
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
					if(LSPINORB)
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

	const double pref = FOUR_PI / sqrt(cell.omega);
	this->tab.zero_out();
	ofs_running<<"\n Init Non-Local PseudoPotential table : ";
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
			for (int iq=0; iq<NQX; iq++)  
			{
				const double q = iq * DQ;
				Mathzone::Spherical_Bessel(kkbeta, cell.atoms[it].r, q, l, jl);

				for (int ir = 0;ir < kkbeta;ir++)
				{
					aux[ir] = cell.atoms[it].betar(ib, ir) *
					          jl[ir] * cell.atoms[it].r[ir];
				} 
				double vqint;
				Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].rab, vqint);
				this->tab(it, ib, iq) = vqint * pref;
			} 
		} 
		delete[] aux;
		delete[] jl;
	}
	timer::tick("ppcell_vnl","init_vnl");
	ofs_running << "\n Init Non-Local-Pseudopotential done." << endl;
	return;
}


complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M)   // pengfei Li  2018-3-23
{
	complex<double> cf;
	if(alpha == 0)
	{
		cf = -sqrt(4*PI/3)*CG(lu,mu,1,1,L,M);
	}
	else if(alpha == 1)
	{
		cf = -sqrt(4*PI/3)*CG(lu,mu,1,2,L,M);
	}
	else if(alpha == 2)
	{
		cf = sqrt(4*PI/3)*CG(lu,mu,1,0,L,M);
	}
	else
	{
		WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
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
	if(test_pp) TITLE("pseudopot_cell_vnl","getvnl_alpha");
	timer::tick("pp_cell_vnl","getvnl_alpha");

	if(lmaxkb < 0) 
	{
		return;
	}
	
	const int npw = kv.ngk[ik];
	int ig, ia, nb, ih, lu, mu;

	double *vq = new double[npw];
	const int x1= (lmaxkb + 2)*(lmaxkb + 2);

	matrix ylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++) 
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}

	vkb1_alpha = new complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		vkb1_alpha[i] = new complex<double>*[nhm];
		for(int j=0; j<nhm; j++)
		{
			vkb1_alpha[i][j] = new complex<double>[npw];
		}
	}	
	
	vkb_alpha = new complex<double>**[3];
	for(int i=0; i<3; i++)
	{
		vkb_alpha[i] = new complex<double>*[nkb];
		for(int j=0; j<nkb; j++)
		{
			vkb_alpha[i][j] = new complex<double>[wf.npwx];
		}
	}
	
	Mathzone::Ylm_Real(x1, npw, gk, ylm);

	MGT.init_Gaunt_CH( lmaxkb + 2 );
	MGT.init_Gaunt( lmaxkb + 2 );
	
	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = ucell.atoms[it].nbeta;
		const int nh = ucell.atoms[it].nh;

		if(test_pp>1) OUT("nbeta",nbeta);

		for(int i=0; i<3; i++)
			for(int j=0; j<nhm; j++)
			{
				ZEROS(vkb1_alpha[i][j], npw);
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
					const double gnorm = gk[ig].norm() * ucell.tpiba;
					vq [ig] = Mathzone::Polynomial_Interpolation(
							this->tab_alpha, it, nb, L, NQX, DQ, gnorm);
					
					for (int M=0; M<2*L+1; M++)
					{
						int lm = L*L + M;
						for (int alpha=0; alpha<3; alpha++)
						{
							complex<double> c = Cal_C(alpha,lu, mu, L, M);
							/*if(alpha == 0)
							{
								cout<<"lu= "<<lu<<"  mu= "<<mu<<"  L= "<<L<<"  M= "<<M<<" alpha = "<<alpha<<"  "<<c<<endl;
							}*/
							vkb1_alpha[alpha][ih][ig] += c * vq[ig] * ylm(lm, ig) * pow( NEG_IMAG_UNIT, L);
						}	
					}
				}
			}
		} // end nbeta

		for (ia=0; ia<ucell.atoms[it].na; ia++) 
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
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
	timer::tick("pp_cell_vnl","getvnl_alpha");
	return;
} 

void pseudopot_cell_vnl::init_vnl_alpha(void)          // pengfei Li 2018-3-23
{
	if(test_pp) TITLE("pseudopot_cell_vnl","init_vnl_alpha");
	timer::tick("ppcell_vnl","init_vnl_alpha");

	for(int it=0;it<ucell.ntype;it++)
	{
		int BetaIndex=0;
		//const int Nprojectors = ucell.atoms[it].nh;
		for (int ib=0; ib<ucell.atoms[it].nbeta; ib++)
		{
			const int l = ucell.atoms[it].lll [ib];
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

	const double pref = FOUR_PI / sqrt(ucell.omega);
	this->tab_alpha.create(ucell.ntype, nbrx, lmaxkb+2, NQX);
	this->tab_alpha.zero_out();
	ofs_running<<"\n Init Non-Local PseudoPotential table( including L index) : ";
	for (int it = 0;it < ucell.ntype;it++)  
	{
		const int nbeta = ucell.atoms[it].nbeta;
		int kkbeta = ucell.atoms[it].kkbeta;

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
				for (int iq = 0; iq < NQX; iq++)
				{
					const double q = iq * DQ;
					Mathzone::Spherical_Bessel(kkbeta, ucell.atoms[it].r, q, L, jl);
					
					for (int ir = 0;ir < kkbeta;ir++)
					{
						aux[ir] = ucell.atoms[it].betar(ib, ir) * jl[ir] * 
								  ucell.atoms[it].r[ir] * ucell.atoms[it].r[ir];
					}
					double vqint;
					Integral::Simpson_Integral(kkbeta, aux, ucell.atoms[it].rab, vqint);
					this->tab_alpha(it, ib, L, iq) = vqint * pref;
				}
			}
		} 
		delete[] aux;
		delete[] jl;
	}
	timer::tick("ppcell_vnl","init_vnl_alpha");
	ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << endl;
	return;
}



void pseudopot_cell_vnl::print_vnl(ofstream &ofs)
{
	out.printr3_d(ofs, " tab : ", tab);
}



int pseudopot_cell_vnl::calculate_nqx(const double &ecutwfc,const double &dq)
{
	int points_of_table = static_cast<int>( sqrt(ecutwfc)/dq + 4 ) ;//* cell_factor;
//	cout<<"\n nqx = "<< points_of_table;
//----------------------------------------------------------
// EXPLAIN : Plus 1 because the formula is transfered from 
// fortran code.
//----------------------------------------------------------
	return points_of_table + 1; 
}
