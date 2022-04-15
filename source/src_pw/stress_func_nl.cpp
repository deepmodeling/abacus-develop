#include "stress_func.h"
#include "../module_base/math_polyint.h"
#include "../module_base/math_ylmreal.h"
#include "../module_base/timer.h"
#include "global.h"

//calculate the nonlocal pseudopotential stress in PW
void Stress_Func::stress_nl(ModuleBase::matrix& sigma)
{
	ModuleBase::TITLE("Stress_Func","stres_nl");
	ModuleBase::timer::tick("Stress_Func","stres_nl");
	
	const int nkb = GlobalC::ppcell.nkb;
	if(nkb == 0) 
	{
		ModuleBase::timer::tick("Stress_Func","stres_nl");
		return;
	}
	double sigmanlc[3][3];
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigmanlc[l][m]=0.0;
		}
	}
	
	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	ModuleBase::ComplexMatrix dbecp( nkb, GlobalV::NBANDS);
	ModuleBase::ComplexMatrix becp( nkb, GlobalV::NBANDS);

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ModuleBase::ComplexMatrix vkb1( nkb, GlobalC::wf.npwx );
	ModuleBase::ComplexMatrix vkb0[3];
	for(int i=0;i<3;i++){
		vkb0[i].create(nkb, GlobalC::wf.npwx);
	}
	ModuleBase::ComplexMatrix vkb2( nkb, GlobalC::wf.npwx );
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {   	  
		if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		GlobalC::wf.npw = GlobalC::kv.ngk[ik];
		// generate vkb
		if (GlobalC::ppcell.nkb > 0)
		{
			GlobalC::ppcell.getvnl(ik);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        ModuleBase::timer::tick("Stress", "cal_becp");
        becp.zero_out();
		for (int ib=0; ib<GlobalV::NBANDS; ib++)
		{
			///
			///only occupied band should be calculated.
			///
			if(GlobalC::wf.wg(ik, ib) < ModuleBase::threshold_wg) continue;
			for (int i = 0; i < nkb; i++) 
			{
                for (int ig = 0; ig < GlobalC::wf.npw; ig++) {
                    becp(i, ib) += GlobalC::wf.evc[ik](ib, ig) *
                                           conj(GlobalC::ppcell.vkb(i, ig));
                }
            }
        }
        ModuleBase::timer::tick("Stress", "cal_becp");
		//becp calculate is over , now we should broadcast this data.
        ModuleBase::timer::tick("Stress", "reduce_complex_double_pool");
		Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);
		ModuleBase::timer::tick("Stress", "reduce_complex_double_pool");

		ModuleBase::timer::tick("Stress", "get_dvnl1");
		for (int i = 0; i < 3; i++) 
		{
			get_dvnl1(vkb0[i], ik, i);
		}
        ModuleBase::timer::tick("Stress", "get_dvnl1");

        ModuleBase::timer::tick("Stress", "get_dvnl2");
        get_dvnl2(vkb2, ik);
		ModuleBase::timer::tick("Stress", "get_dvnl2");

        ModuleBase::Vector3<double> qvec;
        double qvec0[3];

        for (int ipol = 0; ipol < 3; ipol++) 
		{
            for (int jpol = 0; jpol < ipol + 1; jpol++) 
			{
				dbecp.zero_out();
				vkb1.zero_out();
				ModuleBase::timer::tick("Stress", "get_vkb1");
				for (int i = 0; i < nkb; i++) 
				{
					for (int ig = 0; ig < GlobalC::wf.npw; ig++) 
					{
						qvec = GlobalC::wf.get_1qvec_cartesian(ik, ig);
						qvec0[0] = qvec.x;
						qvec0[1] = qvec.y;
						qvec0[2] = qvec.z;

						vkb1(i, ig) += 0.5 * qvec0[ipol] *
											vkb0[jpol](i, ig) +
										0.5 * qvec0[jpol] *
											vkb0[ipol](i, ig);
					} // end ig
					  
				}//end nkb
				ModuleBase::timer::tick("Stress", "get_vkb1");
				ModuleBase::timer::tick("Stress", "dbecp_noevc");
				ModuleBase::ComplexMatrix dbecp_noevc(nkb, GlobalC::wf.npw);
				for (int i = 0; i < nkb; i++) 
				{
					for (int ig = 0; ig < GlobalC::wf.npw;ig++) 
					{
						// first term
						dbecp_noevc(i, ig) -= 2.0 * conj(vkb1(i, ig));
						// second termi
						if (ipol == jpol)
							dbecp_noevc(i, ig) -= conj( GlobalC::ppcell.vkb(i, ig));
						// third term
						qvec =	GlobalC::wf.get_1qvec_cartesian(ik,ig);
						qvec0[0] = qvec.x;
						qvec0[1] = qvec.y;
						qvec0[2] = qvec.z;
						double qm1;
						if(qvec.norm2() > 1e-16) qm1 = 1.0 / qvec.norm(); 
						else qm1 = 0; 
						dbecp_noevc(i,ig)	-= 2.0 * conj(vkb2(i,ig)) * qvec0[ipol] * 
							qvec0[jpol] * qm1 *	GlobalC::ucell.tpiba;
					} // end ig
				}     // end i
				ModuleBase::timer::tick("Stress", "dbecp_noevc");

				ModuleBase::timer::tick("Stress", "get_dbecp");
				for (int ib=0; ib<GlobalV::NBANDS; ib++)
				{
					///
					///only occupied band should be calculated.
					///
					if(GlobalC::wf.wg(ik, ib) < ModuleBase::threshold_wg) continue;
					for (int i=0; i<nkb; i++)
					{
						for (int ig=0; ig<GlobalC::wf.npw; ig++) 
						{
							//first term
							dbecp(i,ib) += GlobalC::wf.evc[ik](ib ,ig) * dbecp_noevc(i, ig);
						}//end ig
					}//end i
				}//end ib
                ModuleBase::timer::tick("Stress", "get_dbecp");

				//              don't need to reduce here, keep
				//              dbecp different in each
				//              processor, and at last sum up
				//              all the forces.
				//              Parallel_Reduce::reduce_complex_double_pool(
				//              dbecp.ptr, dbecp.ndata);

				//              double *cf = new
				//              double[GlobalC::ucell.nat*3];
				//              ModuleBase::GlobalFunc::ZEROS(cf,
				//              GlobalC::ucell.nat);
				ModuleBase::timer::tick("Stress", "get_final_step");
				for (int ib=0; ib<GlobalV::NBANDS; ib++)
				{
					///
					///only occupied band should be calculated.
					///
					if(GlobalC::wf.wg(ik, ib) < ModuleBase::threshold_wg) continue;
					double fac = GlobalC::wf.wg(ik, ib) * 1.0;
					int iat = 0;
					int sum = 0;
					for (int it=0; it<GlobalC::ucell.ntype; it++)
					{
						const int Nprojs = GlobalC::ucell.atoms[it].nh;
						for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
						{
							for (int ip=0; ip<Nprojs; ip++)
							{
								double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip) ;
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
                ModuleBase::timer::tick("Stress","get_final_step");
            }//end jpol
		}//end ipol
	}// end ik

	// sum up forcenl from all processors
	for(int l=0;l<3;l++){
		for(int m=0;m<3;m++){
			if(m>l) 
			{
				sigmanlc[l][m] = sigmanlc[m][l];
			}
			Parallel_Reduce::reduce_double_all( sigmanlc[l][m] ); //qianrui fix a bug for kpar > 1
		}
	}

//        Parallel_Reduce::reduce_double_all(sigmanl.c, sigmanl.nr * sigmanl.nc);
        
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigmanlc[ipol][jpol] *= 1.0 / GlobalC::ucell.omega;
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
	if(ModuleSymmetry::Symmetry::symm_flag)
	{
		GlobalC::symm.stress_symmetry(sigma, GlobalC::ucell);
	}//end symmetry
	
	//  this->print(GlobalV::ofs_running, "nonlocal stress", stresnl);
	ModuleBase::timer::tick("Stress_Func","stres_nl");
	return;
}
 
void Stress_Func::get_dvnl1
(
	ModuleBase::ComplexMatrix &vkb,
	const int ik,
	const int ipol
)
{
	if(GlobalV::test_pp) ModuleBase::TITLE("Stress_Func","get_dvnl1");

	const int lmaxkb = GlobalC::ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = GlobalC::kv.ngk[ik];
	const int nhm = GlobalC::ppcell.nhm;
	int ig, ia, nb, ih;
	ModuleBase::matrix vkb1(nhm, npw);
	vkb1.zero_out();
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix dylm(x1, npw);
	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}
			   
	dylmr2(x1, npw, gk, dylm, ipol);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].nbeta;
		const int nh = GlobalC::ucell.atoms[it].nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;

				//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
				//cout << "\n gk.norm = " << gnorm;

				vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
						GlobalC::ppcell.tab, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );

			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == GlobalC::ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( GlobalC::ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = dylm(lm, ig) * vq [ig];
					}

				}

			} // end ih

		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			std::complex<double> *sk = GlobalC::wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				std::complex<double> pref = pow( ModuleBase::NEG_IMAG_UNIT, GlobalC::ppcell.nhtol(it, ih));      //?
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

void Stress_Func::get_dvnl2(ModuleBase::ComplexMatrix &vkb,
		const int ik)
{
	if(GlobalV::test_pp) ModuleBase::TITLE("Stress","get_dvnl2");
//	ModuleBase::timer::tick("Stress","get_dvnl2");

	const int lmaxkb = GlobalC::ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = GlobalC::kv.ngk[ik];
	const int nhm = GlobalC::ppcell.nhm;
	int ig, ia, nb, ih;
	ModuleBase::matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix ylm(x1, npw);
	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}
	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].nbeta;
		const int nh = GlobalC::ucell.atoms[it].nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
	//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
	//cout << "\n gk.norm = " << gnorm;
				vq [ig] = Polynomial_Interpolation_nl(
						GlobalC::ppcell.tab, it, nb, GlobalV::DQ, gnorm );

			} // enddo

							// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == GlobalC::ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( GlobalC::ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			std::complex<double> *sk = GlobalC::wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				std::complex<double> pref = pow( ModuleBase::NEG_IMAG_UNIT, GlobalC::ppcell.nhtol(it, ih));      //?
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
//	ModuleBase::timer::tick("Stress","get_dvnl2");

	return;
}




double Stress_Func::Polynomial_Interpolation_nl
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
    const double &table_interval,
    const double &x                             // input value
)
{

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


	return y;
}

void Stress_Func::dylmr2 (
	const int nylm,
	const int ngy,
	ModuleBase::Vector3<double> *gk,
	ModuleBase::matrix &dylm,
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

	ModuleBase::matrix ylmaux;
	// dg is the finite increment for numerical derivation:
	// dg = delta |G| = delta * sqrt(gg)
	// dgi= 1 /(delta * sqrt(gg))
	// gx = g +/- dg


	ModuleBase::Vector3<double> *gx = new ModuleBase::Vector3<double> [ngy];
	 

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

	ModuleBase::YlmReal::Ylm_Real(nylm, ngy, gx, dylm);
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

	ModuleBase::YlmReal::Ylm_Real(nylm, ngy, gx, ylmaux);


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
