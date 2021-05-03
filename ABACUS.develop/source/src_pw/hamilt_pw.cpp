#include "tools.h"
#include "global.h"
#include "hamilt_pw.h"
#include "src_global/blas_connector.h"
#include "../src_io/optical.h" // only get judgement to calculate optical matrix or not.
#include "src_pw/myfunc.h"

int Hamilt_PW::moved = 0;

Hamilt_PW::Hamilt_PW()
{
    hpsi = new complex<double>[1];
    spsi = new complex<double>[1];
    GR_index = new int[1];
    Bec = new complex<double>[1];
}

Hamilt_PW::~Hamilt_PW()
{
    delete[] hpsi;
    delete[] spsi;
    delete[] GR_index;
    delete[] Bec;
}


void Hamilt_PW::allocate(
	const int &npwx, 
	const int &npol, 
	const int &nkb, 
	const int &nrxx)
{
    TITLE("Hamilt_PW","allocate");

	assert(npwx > 0);
	assert(npol > 0);
	assert(nkb >=0);
	assert(nrxx > 0);

    delete[] hpsi;
    delete[] spsi;
    delete[] GR_index;
    delete[] Bec;

    this->hpsi = new complex<double> [npwx * npol];
    this->spsi = new complex<double> [npwx * npol];
    this->GR_index = new int[nrxx];
    this->Bec = new complex<double> [nkb];

    ZEROS(this->hpsi, npwx * npol);
    ZEROS(this->spsi, npwx * npol);
    ZEROS(this->GR_index, nrxx);

    return;
}


void Hamilt_PW::init_k(const int ik)
{
    TITLE("Hamilt_PW","init_k");
	
	// mohan add 2010-09-30
	// (1) Which spin to use.
	if(NSPIN==2)
	{
		CURRENT_SPIN = kv.isk[ik];
	}

	// (2) Kinetic energy.
	wf.ekin(ik);

	// (3) Take the local potential.
	for (int ir=0; ir<pw.nrxx; ir++)
	{
		pot.vr_eff1[ir] = pot.vr_eff(CURRENT_SPIN, ir);//mohan add 2007-11-12
	}

	// (4) Calculate nonlocal pseudopotential vkb
	//if (ppcell.nkb > 0 && !LINEAR_SCALING) xiaohui modify 2013-09-02
	if(ppcell.nkb > 0 && (BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")) //xiaohui add 2013-09-02. Attention...
	{
		ppcell.getvnl(ik);
	}

	// (5) The number of wave functions.
	wf.npw = kv.ngk[ik];

	// (6) The index of plane waves.
    for (int ig = 0;ig < wf.npw;ig++)
    {
        this->GR_index[ig] = pw.ig2fftw[ wf.igk(ik, ig) ];
    }
    return;
}


//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
void Hamilt_PW::cinitcgg(
    const int ik,
    const int nstart,
    const int n_band,
    const ComplexMatrix &psi,
    ComplexMatrix &evc,
    double *en)
{
    TITLE("Hamilt_PW","cinitcgg");
    timer::tick("Hamilt_PW","cinitcgg",'G');

	assert(nstart!=0);
	assert(n_band!=0);

    ComplexMatrix hc(nstart, nstart);
    ComplexMatrix sc(nstart, nstart);
    ComplexMatrix hvec(nstart,n_band);
	int dmin,dmax;
	const int npw = kv.ngk[ik];
	if(NSPIN != 4)
	{
		dmin= npw;
		dmax = wf.npwx;
	}
	else {
		dmin = wf.npwx*NPOL;
		dmax = wf.npwx*NPOL;
	}
	//qianrui improve this part 2021-3-14
	complex<double> *aux=new complex<double> [dmax*nstart];
	complex<double> *paux = aux;
	complex<double> *ppsi = psi.c;
	//qianrui replace it
	this->h_psi(psi.c, aux, nstart);
	//for(int m=0;m<nstart;++m)
	//{
	//	this->h_psi(ppsi, paux);
	//	paux += dmax;
	//	ppsi += dmax;
	//}
	char trans1 = 'C';
	char trans2 = 'N';
	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ONE,psi.c,&dmax,aux,&dmax,&ZERO,hc.c,&nstart);
	hc=transpose(hc,false);

	zgemm_(&trans1,&trans2,&nstart,&nstart,&dmin,&ONE,psi.c,&dmax,psi.c,&dmax,&ZERO,sc.c,&nstart);
	sc=transpose(sc,false);
	//After psis are strictly normalized, we should use this part. 
	//for(int m=1;m<nstart;++m)
	//{
	//	sc(m,m) == 1;
	//}

	delete []aux;


	/*//qianrui replace this part
	complex<double> **p = new complex<double>*[nstart];
	for(int i=0; i<nstart; i++)
	{
		if(NPOL==2)
		{
			p[i] = new complex<double>[wf.npwx*NPOL];
		}
		else
		{
			p[i] = new complex<double>[wf.npwx];
		}
	}

	for(int i=0; i<nstart; i++)
	{
		for(int j=0; j<wf.npwx; j++)
		{
			p[i][j] = psi(i,j);
		}
	}
	if(NPOL == 2)
	{
		for(int i=0; i<nstart; i++)
		{
			for(int j=0; j<wf.npwx; j++)
			{
				p[i][j+wf.npwx] = psi(i,j+wf.npwx);
			}
		}
	}

    // Set up the Hamiltonian and Overlap matrices
	complex<double>* hpsi = new complex<double>[dmin];
	complex<double>* spsi = new complex<double>[dmin];

    for (int m=0; m<nstart; m++)
    {
        this->h_1psi(dmin, p[m], hpsi, spsi);

        hc(m,m) = this->just_ddot(dmin, p[m], hpsi);
        sc(m,m) = this->just_ddot(dmin, p[m], spsi);

        if (m+1<nstart)
        {
            for (int j = m + 1;j < nstart;j++)
            {
                hc(j,m) = this->just_ddot(dmin, p[j], hpsi);
                sc(j,m) = this->just_ddot(dmin, p[j], spsi);
                hc(m, j) = conj(hc(j, m));
                sc(m, j) = conj(sc(j, m));
            }
        }
    }
	
	delete[] hpsi;
	delete[] spsi;
	for(int i=0; i<nstart; i++)
	{
		delete[] p[i];
	}
	delete[] p;*/

	// Peize Lin add 2019-03-09
	if("lcao_in_pw"==BASIS_TYPE)
	{
		auto add_Hexx = [&](const double alpha)
		{
			for (int m=0; m<nstart; ++m)
			{
				for (int n=0; n<nstart; ++n)
				{
					hc(m,n) += alpha * exx_lip.get_exx_matrix()[ik][m][n];
				}
			}
		};
		if( 5==xcf.iexch_now && 0==xcf.igcx_now )				// HF
		{
			add_Hexx(1);
		}
		else if( 6==xcf.iexch_now && 8==xcf.igcx_now )			// PBE0
		{
			add_Hexx(exx_global.info.hybrid_alpha);
		}
		else if( 9==xcf.iexch_now && 12==xcf.igcx_now )			// HSE
		{
			add_Hexx(exx_global.info.hybrid_alpha);		
		}
	}

	if(NPROC_IN_POOL>1)
	{
		Parallel_Reduce::reduce_complex_double_pool( hc.c, nstart*nstart );
		Parallel_Reduce::reduce_complex_double_pool( sc.c, nstart*nstart );
	}

	

    hm.cdiaghg(nstart, n_band, hc, sc, nstart, en, hvec);


	// Peize Lin add 2019-03-09
	if("lcao_in_pw"==BASIS_TYPE)
	{
		switch(exx_global.info.hybrid_type)
		{
			case Exx_Global::Hybrid_Type::HF:
			case Exx_Global::Hybrid_Type::PBE0:
			case Exx_Global::Hybrid_Type::HSE:
				exx_lip.k_pack->hvec_array[ik] = hvec;
				break;
		}
	}
		
    //=======================
    //diagonize the H-matrix
    //=======================

// for tests
/*
		cout << setprecision(3);
		out.printV3(ofs_running,kv.kvec_c[ik]);
		out.printcm_norm("sc",sc,1.0e-4);
		out.printcm_norm("hvec",hvec,1.0e-4);
		out.printcm_norm("hc",hc,1.0e-4);
		cout << endl;
*/

	cout << setprecision(5);

//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------
/*
	cout << "  hc matrix" << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			double a = hc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			cout << setw(6) << a;
		}
		cout << endl;
	}

	cout << "  sc matrix" << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			double a = sc(i,j).real();
			if(abs(a) < 1.0e-5) a = 0;
			cout << setw(6) << a;
		}
		cout << endl;
	}

	cout << "\n Band Energy" << endl;
	for(int i=0; i<NBANDS; i++)
	{
		cout << " e[" << i+1 << "]=" << en[i] * Ry_to_eV << endl;
	}
*/
//--------------------------
// KEEP THIS BLOCK FOR TESTS
//--------------------------


	if((BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") && CALCULATION=="nscf" && !Optical::opt_epsilon2)
	{
		ofs_running << " Not do zgemm to get evc." << endl;
	}
	else if((BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") 
		&& ( CALCULATION == "scf" || CALCULATION == "md" || CALCULATION == "relax")) //pengfei 2014-10-13
	{
		// because psi and evc are different here,
		// I think if psi and evc is the same, 
		// there may be problems, mohan 2011-01-01
		char transa = 'N';
		char transb = 'T';
		zgemm_( &transa,
				&transb,
				&dmax, // m: row of A,C
				&n_band, // n: col of B,C
				&nstart, // k: col of A, row of B
				&ONE, // alpha
				psi.c, // A
				&dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
				hvec.c, // B
				&n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
				&ZERO,  // belta
				evc.c, // C
				&dmax ); // LDC: if(N) max(1, m)
	}
	else
	{
		// As the evc and psi may refer to the same matrix, we first
		// create a temporary matrix to story the result. (by wangjp)
		//qianrui improve this part 2021-3-13
		char transa = 'N';
		char transb = 'T';
		ComplexMatrix evctmp(n_band, dmin,false);
		zgemm_(&transa,&transb,&dmin,&n_band,&nstart,&ONE,psi.c,&dmax,hvec.c,&n_band,&ZERO,evctmp.c,&dmin);
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<dmin; ig++)
			{
				evc(ib,ig) = evctmp(ib,ig);
			}
		}
		
		/*qianrui replace this part 
		ComplexMatrix evctmp(n_band, dmax);
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<wf.npw; ig++)
			{
				complex<double> sum = ZERO;
				complex<double> sum1 = ZERO;
				for(int ib2=0; ib2<nstart; ib2++)
				{
					sum += psi(ib2, ig) * hvec(ib2, ib);
					if(NPOL == 2) 
					{
						sum1 +=  psi(ib2,ig + wf.npwx) * hvec(ib2, ib);
					}
				}
				evctmp(ib,ig) = sum;
				if(NPOL == 2)
				{
					evctmp(ib,ig + wf.npwx) = sum1;
				}
			}
		}
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<wf.npw; ig++)
			{
				evc(ib,ig) = evctmp(ib,ig);
				if(NPOL==2) evc(ib, ig + wf.npwx) = evctmp(ib, ig+wf.npwx);
			}
		}*/

		/*
		for(int ib=0; ib<n_band; ib++)
		{
			for(int ig=0; ig<wf.npw; ig++)
			{
				complex<double> diff = evctmp(ib, ig) - evc(ib, ig);
				if( norm(diff) > 1.0e-5 )
				{
					cout << " npw = " << wf.npw << endl;
					cout << " nstart = " << nstart << endl;
					cout << " nbands = " << n_band << endl;
					cout << " wf.npwx = " << wf.npwx << endl;
					cout << " ib=" << ib << " ig=" << ig << " diff=" << diff << endl;
				}
			}
		}
		*/

	}
    //out.printr1_d("en",en,n_band);

//	cout << "\n bands" << endl;
//	for(int ib=0; ib<n_band; ib++)
//	{
//		cout << " ib=" << ib << " " << en[ib] * Ry_to_eV << endl; 
//	}

    //out.printcm_norm("hvec",hvec,1.0e-8);

    timer::tick("Hamilt_PW","cinitcgg",'G');
    return;
}


void Hamilt_PW::h_1psi( const int npw_in, const complex < double> *psi,
                        complex<double> *hpsi, complex < double> *spsi)
{
    this->h_psi(psi, hpsi);

    for (int i=0;i<npw_in;i++)
    {
        spsi[i] = psi[i];
    }
    return;
}


void Hamilt_PW::s_1psi
(
    const int dim,
    const complex<double> *psi,
    complex<double> *spsi
)
{
    for (int i=0; i<dim; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}


void Hamilt_PW::h_psi(const complex<double> *psi_in, complex<double> *hpsi, const int m)
{
    timer::tick("Hamilt_PW","h_psi",'H');
    int i = 0;
    int j = 0;
    int ig= 0;

	//if(NSPIN!=4) ZEROS(hpsi, wf.npw);
	//else ZEROS(hpsi, wf.npwx * NPOL);//added by zhengdy-soc
	int dmax = wf.npwx * NPOL;

	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	complex<double> *tmhpsi;
	const complex<double> *tmpsi_in;
 	if(T_IN_H)
	{	
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0 ; ib < m; ++ib)
		{
			for(ig = 0;ig < wf.npw; ++ig)
			{
				tmhpsi[ig] = wf.g2kin[ig] * tmpsi_in[ig];
			}
			if(NSPIN==4){
				for(ig=wf.npw; ig < wf.npwx; ++ig)
				{
					tmhpsi[ig] = 0;
				}
				tmhpsi +=wf.npwx;
				tmpsi_in += wf.npwx;
				for (ig = 0;ig < wf.npw ;++ig)
				{
					tmhpsi[ig] = wf.g2kin[ig] * tmpsi_in[ig];
				}
				for(ig=wf.npw; ig < wf.npwx; ++ig)
				{
					tmhpsi[ig] =0;
				}
			}
			tmhpsi += wf.npwx;
			tmpsi_in += wf.npwx;
		}
	}

	//------------------------------------
	//(2) the local potential.
	//-----------------------------------
	timer::tick("Hamilt_PW","vloc",'H');
	if(VL_IN_H)
	{
		tmhpsi = hpsi;
		tmpsi_in = psi_in;
		for(int ib = 0 ; ib < m; ++ib)
		{
			if(NSPIN!=4){
				ZEROS( UFFT.porter, pw.nrxx);
				UFFT.RoundTrip( tmpsi_in, pot.vr_eff1, GR_index, UFFT.porter );
				for (j = 0;j < wf.npw;j++)
				{
					tmhpsi[j] += UFFT.porter[ GR_index[j] ];
				}
			}
			else
			{
				complex<double>* porter1 = new complex<double>[pw.nrxx];
				ZEROS( UFFT.porter, pw.nrxx);
				ZEROS( porter1, pw.nrxx);
				for (int ig=0; ig< wf.npw; ig++)
				{
					UFFT.porter[ GR_index[ig]  ] = tmpsi_in[ig];
					porter1[ GR_index[ig]  ] = tmpsi_in[ig + wf.npwx];
				}
				// (2) fft to real space and doing things.
				pw.FFT_wfc.FFT3D( UFFT.porter, 1);
				pw.FFT_wfc.FFT3D( porter1, 1);
				complex<double> sup,sdown;
				for (int ir=0; ir< pw.nrxx; ir++)
				{
					sup = UFFT.porter[ir] * (pot.vr_eff(0,ir) + pot.vr_eff(3,ir)) +
						porter1[ir] * (pot.vr_eff(1,ir) - complex<double>(0.0,1.0) * pot.vr_eff(2,ir));
					sdown = porter1[ir] * (pot.vr_eff(0,ir) - pot.vr_eff(3,ir)) +
					UFFT.porter[ir] * (pot.vr_eff(1,ir) + complex<double>(0.0,1.0) * pot.vr_eff(2,ir));
					UFFT.porter[ir] = sup;
					porter1[ir] = sdown;
				}
				// (3) fft back to G space.
				pw.FFT_wfc.FFT3D( UFFT.porter, -1);
				pw.FFT_wfc.FFT3D( porter1, -1);

				for (j = 0;j < wf.npw;j++)
				{
					tmhpsi[j] += UFFT.porter[ GR_index[j] ];
				}
				for (j = 0;j < wf.npw;j++ )
				{
					tmhpsi[j+wf.npwx] += porter1[ GR_index[j] ];
				}
				delete[] porter1;
			}
			tmhpsi += dmax;
			tmpsi_in += dmax;
		}
	}
	timer::tick("Hamilt_PW","vloc",'H');

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	timer::tick("Hamilt_PW","vnl",'H');
	if(VNL_IN_H)
	{
		if ( ppcell.nkb > 0)
		{
			//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//qianrui optimize 2021-3-31
			int nkb=ppcell.nkb;
			ComplexMatrix becp(NPOL * m, nkb, false);
			char transa = 'C';
			char transb = 'N';
			if(m==1 && NPOL==1)
			{
				int inc = 1;
				zgemv_(&transa, &wf.npw, &nkb, &ONE, ppcell.vkb.c, &wf.npwx, psi_in, &inc, &ZERO, becp.c, &inc);
			}
			else
			{
				int npm = NPOL * m;
				zgemm_(&transa,&transb,&nkb,&npm,&wf.npw,&ONE,ppcell.vkb.c,&wf.npwx,psi_in,&wf.npwx,&ZERO,becp.c,&nkb);
				//add_vuspsi is moddified, thus tranpose not needed here.
				//if(NONCOLIN)
				//{
				//	ComplexMatrix partbecp(NPOL, nkb ,false);
				//	for(int ib = 0; ib < m; ++ib)
				//	{
//
				//		for ( i = 0;i < NPOL;i++)
				//			for (j = 0;j < nkb;j++)
				//				partbecp(i, j) = tmbecp[i*nkb+j];
				//		for (j = 0; j < nkb; j++)
				//			for (i = 0;i < NPOL;i++)
				//				tmbecp[j*NPOL+i] = partbecp(i, j);
				//		tmbecp += NPOL * nkb;
				//	}
				//}
			}

			Parallel_Reduce::reduce_complex_double_pool( becp.c, nkb * NPOL * m);
			this->add_vuspsi(hpsi, becp.c, m);
			//======================================================================
			/*complex<double> *becp = new complex<double>[ ppcell.nkb * NPOL ];
			ZEROS(becp,ppcell.nkb * NPOL);
			for (i=0;i< ppcell.nkb;i++)
			{
				const complex<double>* p = &ppcell.vkb(i,0);
				const complex<double>* const p_end = p + wf.npw;
				const complex<double>* psip = psi_in; 
				for (;p<p_end;++p,++psip)
				{
					if(!NONCOLIN) becp[i] += psip[0]* conj( p[0] );
					else{
						becp[i*2] += psip[0]* conj( p[0] );
						becp[i*2+1] += psip[wf.npwx]* conj( p[0] );
					}
				}
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, ppcell.nkb * NPOL);
			this->add_vuspsi(hpsi, becp);
			delete[] becp;*/
			//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		}
	}
	timer::tick("Hamilt_PW","vnl",'H');
    timer::tick("Hamilt_PW","h_psi",'H');
    return;
}

void Hamilt_PW::add_vuspsi(complex<double> *hpsi_in,const complex<double> *becp, const int m)
{
    timer::tick("Hamilt_PW","add_vuspsi",'I');
	int nkb = ppcell.nkb;
	complex<double> *Ps  = new complex<double> [nkb * NPOL * m];
    ZEROS( Ps, NPOL * m * nkb);

    int sum = 0;
    int iat = 0;
    // this function sum up each non-local pseudopotential located in each atom,
    // all we need to do is put the right Dij coefficient to each becp, which
    // is calculated before.
    if(NSPIN!=4)
	{
		for (int it=0; it<ucell.ntype; it++)
		{
			const int Nprojs = ucell.atoms[it].nh;
			for (int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				// each atom has Nprojs, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<Nprojs; ip++)
				{
					for (int ip2=0; ip2<Nprojs; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							Ps[(sum + ip2) * m + ib] += ppcell.deeq(CURRENT_SPIN, iat, ip, ip2) * becp[ib * nkb + sum + ip];
						}//end ib
					}// end ih
				}//end jh 
				sum += Nprojs;
				++iat;
			} //end na
		} //end nt
	}
	else
	{
		for (int it=0; it<ucell.ntype; it++)
		{
			int psind,becpind;
			complex<double> becp1,becp2;
			const int Nprojs = ucell.atoms[it].nh;
			for (int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				// each atom has Nprojs, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<Nprojs; ip++)
				{
					for (int ip2=0; ip2<Nprojs; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							psind = (sum+ip2) * 2 * m + ib * 2;
							becpind = ib*nkb*2 + sum + ip;
							becp1 =  becp[becpind];
							becp2 =  becp[becpind + nkb];
							Ps[psind] += ppcell.deeq_nc(0, iat, ip2, ip) * becp1
								+ppcell.deeq_nc(1, iat, ip2, ip) * becp2;
							Ps[psind +1] += ppcell.deeq_nc(2, iat, ip2, ip) * becp1
								+ppcell.deeq_nc(3, iat, ip2, ip) * becp2;
						}//end ib
					}// end ih
				}//end jh
				sum += Nprojs;
				++iat;
			} //end na
		} //end nt
	}

	/*
    for (int ig=0;ig<wf.npw;ig++)
    {
        for (int i=0;i< ppcell.nkb;i++)
        {
            hpsi_in[ig]+=Ps[i]*ppcell.vkb(i,ig);
        }
    }
	*/


	// use simple method.
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//qianrui optimize 2021-3-31
	char transa = 'N';
	char transb = 'T';
	if(NPOL==1 && m==1)
	{
		int inc = 1;
		zgemv_(&transa, &wf.npw, &ppcell.nkb, &ONE, ppcell.vkb.c, &wf.npwx, Ps, &inc, &ONE, hpsi_in, &inc);
	}
	else
	{
		int npm = NPOL*m;
		zgemm_(&transa,&transb,&wf.npw,&npm,&ppcell.nkb,&ONE,ppcell.vkb.c,&wf.npwx,Ps,&npm,&ONE,hpsi_in,&wf.npwx);
	}
	//======================================================================
	/*if(!NONCOLIN)
	for(int i=0; i<ppcell.nkb; i++)
	{
		complex<double>* p = &ppcell.vkb(i,0);
		complex<double>* p_end = p + wf.npw;
		complex<double>* hp = hpsi_in;
		complex<double>* psp = &Ps[i];
		for (;p<p_end;++p,++hp)
		{
			hp[0] += psp[0] * p[0];
		}
	}
	else
	for(int i=0; i<ppcell.nkb; i++)
	{
		complex<double>* p = &ppcell.vkb(i,0);
		complex<double>* p_end = p + wf.npw;
		complex<double>* hp = hpsi_in;
		complex<double>* hp1 = hpsi_in + wf.npwx;
		complex<double>* psp = &Ps[i*2];
		for (;p<p_end;p++,++hp,++hp1)
		{
			hp[0] += psp[0] * (p[0]);
			hp1[0] += psp[1] * (p[0]);
		}
	}*/
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	delete []Ps;
    timer::tick("Hamilt_PW","add_vuspsi",'I');
    return;
}

void Hamilt_PW::diag_zheev(const int &npw_in, ComplexMatrix &psi, const int &nband, double *em, double *err)
{
    TITLE("Hamilt_PW","diag_zheev");
    assert(nband < npw_in) ;

    // if flag =0, this means success.
    // RedM means Reduced matrix, because the dimension of plane wave
    // is much larger than nbands.

    ComplexMatrix RedM(nband, nband);
    complex<double> * eta =  new complex<double>[npw_in] ;
    complex<double> * hpsi1 =  new complex<double>[npw_in] ;
    complex<double> * spsi1 =  new complex<double>[npw_in] ;

	ZEROS(eta, npw_in);
	ZEROS(hpsi1, npw_in);
	ZEROS(spsi1, npw_in);

    double * tmpen = new double[nband] ;

    assert(eta  !=  0) ;
    assert(hpsi1  !=    0) ;
    assert(spsi1    !=  0) ;
    assert(tmpen != 0) ;

    // <j|H|i>, where j < = i is calculated.
    // first calculate eta =|i>, and hpsi = H|i>

    complex<double> tmp ;
    complex<double> tmp1 ;

    // calculate tmpen[i]

    for (int i = 0; i < nband ; i++)
    {
        dcopy(psi, i, eta);
        h_1psi(npw_in, eta, hpsi1, spsi1) ;

        tmp = ZERO ;
        tmp1 = ZERO ;

        for (int ig = 0; ig < npw_in; ig++)
        {
            tmp += conj(eta[ig]) * hpsi1[ig] ;
            tmp1 += conj(eta[ig]) * eta[ig] ;
        }

        tmp /= tmp1 ;

        tmpen[i] = tmp.real() ;

        for (int j = 0; j <= i; j++)
        {
            // calculate H[i,j] = <i|eta> where |eta> = H|j>
            // RedM(j, i) = <j|H|i>

            tmp = ZERO ;

            for (int ig = 0; ig < npw_in; ig++)
            {
                tmp += conj(psi(j, ig)) * hpsi1[ig];
            }

            RedM(j, i) = tmp ;

            if (i != j) RedM(i, j) = conj(tmp) ;
        }
    }

    // This is for calling zheev.
    char  jobz = 'V' ;  // eigenvalues and eigenvectors

    char  uplo = 'U' ;  // the upper is stored.

    const int lwork = 2 * nband;

    complex<double> * work = new complex<double>[lwork]() ;

    double * rwork = new double[3*nband-2] ;

    int info = 0 ;

    // diag by calling zheev
    // basis are psi_0, psi_1, psi_2...

    LapackConnector::zheev(jobz, uplo, nband, RedM, nband, em, work, lwork, rwork, &info);

    delete[] eta  ;

    delete[] hpsi1 ;

    delete[] spsi1 ;

    //  change back to plane wave basis.

//  cout << " calling zheev is performed " << endl ;
//  cout << " zheev output info... " << endl;

    // delete the allocated data array.

    delete[] work ;

    delete[] rwork ;

    // cout the infomation from zheev...

    // cout the infomation from zheev...
    if (info == 0)
    {
        cout << "  successful exit of zheev " << endl ;
    }
    else if (info < 0)
    {
        cout << " the i-th argument had an illegal value. info =  " << info << endl ;
    }
    else
    {
        cout << "the algorithm failed to converge. info = " << info << endl ;
    }

    ComplexMatrix kpsi(nband, npw_in);

    // kpsi = c_1 \psi_1 + c_2 \psi_2 + \cdots + c_N \psi_N
    // For the mth wavefunction, |psi(m)> = U(1, m) \psi_1 + ... U(k, m) \psi_k + ...
    // So, |psi(m)>_j means the jth element of |psi(m)>
    //  We have |psi(m)>_j =  U(1, m) \psi(1, j) + U(2, m) \psi(2, j) + ...
    // = \sum_{k=1}^{nband} U(k, m) \psi(k, j)

    // Store the wavefunction in kpsi

    for (int m = 0; m < nband; m++)
    {
        for (int j = 0; j < npw_in; j++)
        {
            tmp = ZERO ;

            for (int  k = 0; k < nband; k++)
            {
                tmp += RedM(k, m) * psi(k, j);
            }

            kpsi(m, j) = tmp ;
        }
    }


    // update the wavefunction of psi
    for (int m = 0; m < nband ; m++)
    {
        for (int ig = 0; ig < npw_in; ig++)
        {
            psi(m, ig) = kpsi(m, ig);
        }
    }

    // calculate error of the last results.
    // err = || H\psi -E\psi||

    cout << " callilng cal_err " << endl ;

    cal_err(npw_in, psi, nband, em, err);


    // out put the results.
    cout<<setw(6)<<"Bands"
        <<setw(12)<<"energy(ev)"
        <<setw(12)<<"err"
        <<setw(25)<<"||H * psi - E * psi ||\n";

    for (int m = 0; m < 5; m++)
    {
        cout << setw(6) << m
             << setw(12) << em[m] * Ry_to_eV
             << setw(12) << err[m]
             << setw(25) << tmpen[m] * Ry_to_eV << endl ;
    }

    cout << " end of diag_zheev " << endl ;

    return;
}

void Hamilt_PW::cal_err
(
    const int &npw_in,
    ComplexMatrix &psi,
    const int &nband,
    double *em,
    double *err
)
{
//	TITLE("Hamilt_PW", "cal_err");
//	cout << "\n npw_in = " << npw_in << endl;

    assert(nband < npw_in);
    timer::tick("Hamilt_PW", "cal_err") ;

    complex<double> *psitmp =  new complex<double>[npw_in]();
    complex<double> *hpsitmp =  new complex<double>[npw_in]();
    complex<double> *spsitmp =  new complex<double>[npw_in]();

    complex<double> tmp1 ;

    for (int m = 0; m < nband; m++)
    {
//		cout << "\n m = " << m << endl;
        dcopy(psi, m, psitmp) ;
        h_1psi(npw_in, psitmp, hpsitmp, spsitmp);

        complex<double> tmp = ZERO;

        for (int ig=0;  ig<npw_in; ig++)
        {
            tmp1 =  hpsitmp[ig] - em[m] * psitmp[ig];
            tmp +=  conj(tmp1) * tmp1;
        }

        // err[m] = ||H\psitmp - \lambda_m \psitmp||
        err[m] = sqrt( tmp.real() ) ;
    }

    delete[] psitmp;
    delete[] hpsitmp;
    delete[] spsitmp;

//	cout << " calculate error of the wavefunctions " << endl ;
    timer::tick("Hamilt_PW", "cal_err") ;
    return;
}

double Hamilt_PW::ddot_real
(
    const int &dim,
    const complex<double>* psi_L,
    const complex<double>* psi_R
)const
{
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //qianrui modify 2021-3-14
    //Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
    int dim2=2*dim;
    double *pL,*pR;
    pL=(double *)psi_L;
    pR=(double *)psi_R;
    double result=LapackConnector::dot(dim2,pL,1,pR,1);
    Parallel_Reduce::reduce_double_pool( result );
    return result;
    //======================================================================
    /*complex<double> result(0,0);
    for (int i=0;i<dim;i++)
    {
        result += conj( psi_L[i] ) * psi_R[i];
    }
    Parallel_Reduce::reduce_complex_double_pool( result );
    return result.real();*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
}

complex<double> Hamilt_PW::ddot(
    const int & dim,
    const complex<double> * psi_L,
    const complex<double> * psi_R
)const
{
	complex<double> result = ZERO;
	const int incx = 1;
	const int incy = 1;
	// mohan add 2010-10-11
	zdotc_(&result, &dim, psi_L, &incx, psi_R, &incy); 
    
	if(NPROC_IN_POOL>1)
	{
		Parallel_Reduce::reduce_complex_double_pool( result );
	}
    return result;
}

complex<double> Hamilt_PW::just_ddot(
    const int & dim,
    const complex<double> * psi_L,
    const complex<double> * psi_R
)const
{
	complex<double> result = ZERO;

	// mohan add 2010-10-11
//	zdotc_(&result, &dim, psi_L, &incx, psi_R, &incy);  

	// mohan update 2011-09-21
	static int warn_about_zdotc=true;
	if(warn_about_zdotc)
	{
		ofs_warning << " in Hamilt_PW::just_ddot, sometimes zdotc is not available due to GNU compiler!!!" << endl;
		ofs_warning << " So here I use simple for cicle to replace zdotc, but it will affect the speed." << endl;
		warn_about_zdotc=false;
	}
	for(int i=0; i<dim; ++i)
	{
		result += conj(psi_L[i])*psi_R[i];	
	}

    return result;
}



// this return <psi(m)|psik>
complex<double> Hamilt_PW::ddot(
    const int & dim,
    const ComplexMatrix &psi,
    const int & m,
    const complex<double> *psik
)const
{
    complex<double> result(0, 0);
    assert(dim > 0) ;

    for (int i = 0; i < dim ; i++)
    {
        result += conj(psi(m, i)) *  psik[i] ;
    }

    Parallel_Reduce::reduce_complex_double_pool( result );

    return result;
}  // end of ddot

