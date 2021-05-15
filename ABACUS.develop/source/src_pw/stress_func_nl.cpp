#include "stress_func.h"
#include "../src_global/math_polyint.h"
#include "../src_global/math_ylmreal.h"

//calculate the nonlocal pseudopotential stress in PW
void Stress_Func::stress_nl(matrix& sigma){
	TITLE("Stress_Func","stres_nl");
	timer::tick("Stress_Func","stres_nl",'F');
	
	const int nkb = ppcell.nkb;
	if(nkb == 0) return;

	double sigmanlc[3][3];
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigmanlc[l][m]=0.0;
		}
	}
	
	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	ComplexMatrix dbecp( nkb, NBANDS);
	ComplexMatrix becp( nkb, NBANDS);

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ComplexMatrix vkb1( nkb, wf.npwx );
	ComplexMatrix vkb0[3];
	for(int i=0;i<3;i++){
		vkb0[i].create(nkb, wf.npwx);
	}
	ComplexMatrix vkb2( nkb, wf.npwx );
    for (int ik = 0;ik < kv.nks;ik++)
    {
		for(int i=0;i<3;i++){
			vkb0[i].zero_out();
		}
		vkb2.zero_out();      
		  
		if (NSPIN==2) CURRENT_SPIN = kv.isk[ik];
		wf.npw = kv.ngk[ik];
		// generate vkb
		if (ppcell.nkb > 0)
		{
			ppcell.getvnl(ik);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
		becp.zero_out();
		for (int ib=0; ib<NBANDS; ib++)
		{
			for (int i=0;i<nkb;i++)
			{
				for (int ig=0; ig<wf.npw; ig++)
				{
					becp(i,ib) += wf.evc[ik](ib,ig)* conj( ppcell.vkb(i,ig) );
				}
			}
		}
		Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);
			
		for(int i=0;i<3;i++) {
			get_dvnl1(vkb0[i],ik,i);
		}
				
		get_dvnl2(vkb2,ik);

		Vector3<double> qvec;
		double qvec0[3];
				
		for (int ipol = 0; ipol<3; ipol++)
		{
			for(int jpol = 0; jpol < ipol+1; jpol++)
			{
				dbecp.zero_out();    
				vkb1.zero_out();
				for (int i = 0;i < nkb;i++)
				{
					for (int ig=0; ig<wf.npw; ig++)  
					{
						qvec = wf.get_1qvec_cartesian(ik,ig) ;
						qvec0[0] = qvec.x;
						qvec0[1] = qvec.y;
						qvec0[2] = qvec.z;
					  
						vkb1(i, ig) += 0.5 * qvec0[ipol] * vkb0[jpol](i,ig)
								   + 0.5 * qvec0[jpol] * vkb0[ipol](i,ig) ;
					}//end ig
							
							
				/*  if (kpol==0)
					{
						for (int ig=0; ig<wf.npw; ig++)
						{
							qvec = wf.get_1qvec_cartesian(ik,ig);
							vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.x;
						}
					}
					if (kpol==1)
					{
						for (int ig=0; ig<wf.npw; ig++)
						{
							qvec = wf.get_1qvec_cartesian(ik,ig);
							vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.y;
						}
					}
					if (kpol==2)
					{
						for (int ig=0; ig<wf.npw; ig++)
						{
							qvec = wf.get_1qvec_cartesian(ik,ig);
							vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * qvec.z;
						}
					}*/
					  
				}//end nkb
 

				for (int ib=0; ib<NBANDS; ib++)
				{
					for (int i=0; i<nkb; i++)
					{
						for (int ig=0; ig<wf.npw; ig++)
						{
							//first term
							dbecp(i,ib) = dbecp(i,ib) - 2.0 * wf.evc[ik](ib,ig) * conj( vkb1(i,ig) ) ;
							//second termi
							if(ipol == jpol)
								 dbecp(i,ib) += -1.0 * wf.evc[ik](ib,ig)* conj( ppcell.vkb(i,ig) );
							//third term
							qvec = wf.get_1qvec_cartesian(ik,ig);
							qvec0[0] = qvec.x;
							qvec0[1] = qvec.y;
							qvec0[2] = qvec.z;
							double qm1; 
							if(qvec.norm() > 1e-8) qm1 = 1.0 / qvec.norm();
							else qm1 = 0;
							dbecp(i,ib) +=  -2.0 * wf.evc[ik](ib,ig) * conj(vkb2(i,ig)) * qvec0[ipol] * qvec0[jpol] * qm1 * ucell.tpiba;
						}//end ig
					}//end i
				}//end ib

//              don't need to reduce here, keep dbecp different in each processor,
//              and at last sum up all the forces.
//              Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//              double *cf = new double[ucell.nat*3];
//              ZEROS(cf, ucell.nat);
				for (int ib=0; ib<NBANDS; ib++)
				{
					double fac = wf.wg(ik, ib) * 1.0;
					int iat = 0;
					int sum = 0;
					for (int it=0; it<ucell.ntype; it++)
					{
						const int Nprojs = ucell.atoms[it].nh;
						for (int ia=0; ia<ucell.atoms[it].na; ia++)
						{
							for (int ip=0; ip<Nprojs; ip++)
							{
								double ps = ppcell.deeq(CURRENT_SPIN, iat, ip, ip) ;
								const int inkb = sum + ip;
								//out<<"\n ps = "<<ps;

							 
								const double dbb = ( conj( dbecp( inkb, ib) ) * becp( inkb, ib) ).real();
								sigmanlc[ipol][ jpol] -= ps * fac * dbb;
							 
							}//end ip
							++iat;        
							sum+=Nprojs;
						}//ia
					} //end it
				} //end band
			}//end jpol
		}//end ipol
	}// end ik

	// sum up forcenl from all processors
	for(int l=0;l<3;l++){
		for(int m=0;m<3;m++){
			if(m>l) sigmanlc[l][m] = sigmanlc[m][l];
				Parallel_Reduce::reduce_double_pool( sigmanlc[l][m] );
		}
	}

//        Parallel_Reduce::reduce_double_all(sigmanl.c, sigmanl.nr * sigmanl.nc);
        
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigmanlc[ipol][jpol] *= 1.0 / ucell.omega;
		}
	}
	
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigma(ipol,jpol) = sigmanlc[ipol][jpol] ;
		}
	}
	//do symmetry
	if(Symmetry::symm_flag)
	{
		symm.stress_symmetry(sigma);
	}//end symmetry
	
	//  this->print(ofs_running, "nonlocal stress", stresnl);
	timer::tick("Stress_Func","stres_nl");
	return;
}
 
void Stress_Func::get_dvnl1
(
	ComplexMatrix &vkb,
	const int ik,
	const int ipol
)
{
	if(test_pp) TITLE("Stress_Func","get_dvnl1");

	const int lmaxkb = ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = kv.ngk[ik];
	const int nhm = ppcell.nhm;
	int ig, ia, nb, ih;
	matrix vkb1(nhm, npw);
	vkb1.zero_out();
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	matrix dylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}
			   
	dylmr2(x1, npw, gk, dylm, ipol);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
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

				//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
				//cout << "\n gk.norm = " << gnorm;

				vq [ig] = PolyInt::Polynomial_Interpolation(
						ppcell.tab, it, nb, NQX, DQ, gnorm );

			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = dylm(lm, ig) * vq [ig];
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
				complex<double> pref = pow( NEG_IMAG_UNIT, ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
		delete [] sk;
		} // end ia
	} // enddo
	delete [] gk;
	delete [] vq;
	return;
}//end get_dvnl1

void Stress_Func::get_dvnl2(ComplexMatrix &vkb,
		const int ik)
{
	if(test_pp) TITLE("Stress","get_dvnl2");
//	timer::tick("Stress","get_dvnl2");

	const int lmaxkb = ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = kv.ngk[ik];
	const int nhm = ppcell.nhm;
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
	YlmReal::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
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
	//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
	//cout << "\n gk.norm = " << gnorm;
				vq [ig] = Polynomial_Interpolation_nl(
						ppcell.tab, it, nb, DQ, gnorm );

			} // enddo

							// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( ppcell.nhtolm(it, ih) );
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
				complex<double> pref = pow( NEG_IMAG_UNIT, ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
			++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	}	 // enddo

	delete [] gk;
	delete [] vq;
//	timer::tick("Stress","get_dvnl2");

	return;
}




double Stress_Func::Polynomial_Interpolation_nl
(
    const realArray &table,
    const int &dim1,
    const int &dim2,
    const double &table_interval,
    const double &x                             // input value
)
{
//      timer::tick("Mathzone","Poly_Interpo_2");
	assert(table_interval>0.0);
	const double position = x  / table_interval;
	const int iq = static_cast<int>(position);

	const double x0 = position - static_cast<double>(iq);
	const double x1 = 1.0 - x0;
	const double x2 = 2.0 - x0;
	const double x3 = 3.0 - x0;
	const double y=
			( table(dim1, dim2, iq)   * (-x2*x3-x1*x3-x1*x2) / 6.0 +
			table(dim1, dim2, iq+1) * (+x2*x3-x0*x3-x0*x2) / 2.0 -
			table(dim1, dim2, iq+2) * (+x1*x3-x0*x3-x0*x1) / 2.0 +
			table(dim1, dim2, iq+3) * (+x1*x2-x0*x2-x0*x1) / 6.0 )/table_interval ;

//      timer::tick("Mathzone","Poly_Interpo_2");
	return y;
}

void Stress_Func::dylmr2 (
	const int nylm,
	const int ngy,
	Vector3<double> *gk,
	matrix &dylm,
	const int ipol)
{
  //-----------------------------------------------------------------------
  //
  //     compute \partial Y_lm(G) \over \partial (G)_ipol
  //     using simple numerical derivation (SdG)
  //     The spherical harmonics are calculated in ylmr2
  //
  //int nylm, ngy, ipol;
  // number of spherical harmonics
  // the number of g vectors to compute
  // desired polarization
  //double g (3, ngy), gg (ngy), dylm (ngy, nylm)
  // the coordinates of g vectors
  // the moduli of g vectors
  // the spherical harmonics derivatives
  //
	int ig, lm;
	// counter on g vectors
	// counter on l,m component

	const double delta = 1e-6;
	double *dg, *dgi;

	matrix ylmaux;
	// dg is the finite increment for numerical derivation:
	// dg = delta |G| = delta * sqrt(gg)
	// dgi= 1 /(delta * sqrt(gg))
	// gx = g +/- dg


	Vector3<double> *gx = new Vector3<double> [ngy];
	 

	dg = new double [ngy];
	dgi = new double [ngy];

	ylmaux.create (nylm, ngy);

	dylm.zero_out();
	ylmaux.zero_out();

	for( ig = 0;ig< ngy;ig++){
		gx[ig] = gk[ig];
	}
	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for( ig = 0;ig< ngy;ig++){
		dg [ig] = delta * gx[ig].norm() ;
		if (gx[ig].norm2() > 1e-9) {
			dgi [ig] = 1.0 / dg [ig];
		}
		else{
			dgi [ig] = 0.0;
		}
	}
	//$OMP END PARALLEL DO

	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for( ig = 0;ig< ngy;ig++)
	{
		if(ipol==0)
			gx [ig].x = gk[ ig].x + dg [ig];
		else if(ipol==1)
			gx [ig].y = gk [ ig].y + dg [ig];
		else if(ipol==2)
			gx [ig].z = gk [ ig].z + dg [ig];
	}
	//$OMP END PARALLEL DO

	YlmReal::Ylm_Real(nylm, ngy, gx, dylm);
	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for(ig = 0;ig< ngy;ig++)
	{
		if(ipol==0)
			gx [ig].x = gk [ ig].x - dg [ig];
		else if(ipol==1)
			gx [ig].y = gk [ ig].y - dg [ig];
		else if(ipol==2)
			gx [ig].z = gk [ ig].z - dg [ig];
	}
	//$OMP END PARALLEL DO

	YlmReal::Ylm_Real(nylm, ngy, gx, ylmaux);


	//  zaxpy ( - 1.0, ylmaux, 1, dylm, 1);
	for( lm = 0;lm< nylm;lm++)
	{
		for(ig = 0;ig< ngy;ig++)
		{
			dylm (lm,ig) = dylm(lm,ig) - ylmaux(lm,ig);
		}
	}


	for( lm = 0;lm< nylm;lm++)
	{
		//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
		for(ig = 0;ig< ngy;ig++)
		{
			dylm (lm,ig) = dylm(lm,ig) * 0.5 * dgi [ig];
		}
		//$OMP END PARALLEL DO
	}
	delete[] gx;
	delete[] dg;
	delete[] dgi;

	return;
}
