#include "wf_atomic.h"
#include "global.h"

WF_atomic::WF_atomic()
{
    evc  = new ComplexMatrix[1];
    wanf2= new ComplexMatrix[1];
}

WF_atomic::~WF_atomic()
{ 
	if(test_deconstructor)
	{
		cout << " ~WF_atomic()" << endl;
	}
    delete[] evc;
    delete[] wanf2;
}

//==========================================================
// MEMBER FUNCTION :
// NAME : init_at_1(init a table with the radial Fourier
// transform of the atomic WF_atomictions)
//==========================================================
void WF_atomic::init_at_1(void)
{
    if (test_wf) TITLE("WF_atomic","init_at_1");
    timer::tick("WF_atomic","init_at_1");

	ofs_running << "\n Make real space PAO into reciprocal space." << endl;

    this->print_PAOs();

//----------------------------------------------------------
// EXPLAIN : Find the type of atom that has most mesh points.
//----------------------------------------------------------
    int ndm = 0;
    for (int it=0; it<ucell.ntype; it++)
    {
        ndm = (ucell.atoms[it].msh > ndm) ? ucell.atoms[it].msh : ndm;
    }
	OUT(ofs_running,"max mesh points in Pseudopotential",ndm);

    // needed to normalize atomic wfcs (not a bad idea in general and
    // necessary to compute correctly lda+U projections)
    ppcell.tab_at.zero_out();
//----------------------------------------------------------
// EXPLAIN : If use gauss orbitals to represent aotmic
// orbitals (controlled by parameters)
//
// USE GLOBAL CLASS VARIABLES :
// NAME : ucell.atoms.nchi
// NAME : ucell.atoms.msh(number of mesh points)
// NAME : ucell.atoms.r
//----------------------------------------------------------
    const int startq = 0;
    const double pref = FOUR_PI / sqrt(ucell.omega);
    double *aux = new double[ndm];
    double *vchi = new double[ndm];

	OUT(ofs_running,"dq(describe PAO in reciprocal space)",DQ);
	OUT(ofs_running,"max q",NQX);

    for (int it=0; it<ucell.ntype; it++)
    {
		Atom* atom = &ucell.atoms[it];

		ofs_running <<"\n number of pseudo atomic orbitals for "
		<< atom->label << " is " << atom->nchi << endl;	

        for (int ic=0; ic<atom->nchi ;ic++)
        {
			//cout << "\n T=" << it << " ic=" << ic << endl;

            // check the unit condition
            double *inner_part = new double[atom->msh];
            for (int ir=0; ir<atom->msh; ir++)
            {
                inner_part[ir] = atom->chi(ic,ir) * atom->chi(ic,ir);
            }
            double unit = 0.0;
            Mathzone::Simpson_Integral(atom->msh, inner_part, atom->rab, unit);
            delete[] inner_part;

			ofs_running << " the unit of pseudo atomic orbital is " << unit; 

            //=================================
            // normalize radial wave functions
            //=================================
            for (int ir=0; ir<atom->msh; ir++)
            {
                atom->chi(ic,ir) /= sqrt(unit);
            }

            //===========
            // recheck
            //===========
            inner_part = new double[atom->msh];
            for (int ir=0; ir<atom->msh; ir++)
            {
                inner_part[ir] = atom->chi(ic,ir) * atom->chi(ic,ir);
            }
            unit = 0.0;
            Mathzone::Simpson_Integral(atom->msh, inner_part, atom->rab, unit);
            delete[] inner_part;

			ofs_running << ", renormalize to " << unit << endl;

            if (atom->oc[ic] >= 0.0)
            {
                const int l = atom->lchi[ic];
                for (int iq=startq; iq<NQX; iq++)
                {
                    const double q = DQ * iq;
                    Mathzone::Spherical_Bessel(atom->msh, atom->r, q, l, aux);
                    for (int ir = 0;ir < atom->msh;ir++)
                    {
                        vchi[ir] = atom->chi(ic,ir) * aux[ir] * atom->r[ir];
                    }

                    double vqint = 0.0;
                    Mathzone::Simpson_Integral(atom->msh, vchi, atom->rab, vqint);

                    ppcell.tab_at(it, ic, iq) =  vqint * pref;
                    //				if( it == 0 && ic == 0 )
                    //				{
                    //
                    //					for (ir = 0;ir < ucell.atoms[it].msh;ir++)
                    //						ofs_running << setprecision(20) << "\n vchi(" << ir << ")=" << vchi[ir];
                    //					ofs_running<<"\n aux[0] = "<<aux[0];
                    //					ofs_running<<"\n msh = "<< ucell.atoms[it].msh;
                    //					ofs_running<<"\n tab_at : "<<ppcell.tab_at(it, ic, iq) ;
                    //					ofs_running<<"\n pref = "<<pref;
                    //				}
                } // enddo
            } // endif
        } // enddo
    } // enddo

    delete [] aux;
    delete [] vchi;
    timer::tick("WF_atomic","init_at_1");
    return;
}// end init_at_1

void WF_atomic::print_PAOs(void)const
{
    if (MY_RANK!=0) return;
    for (int it=0; it<ucell.ntype; it++)
    {
        for (int icc=0; icc<ucell.atoms[it].nchi ;icc++)
        {
            int ic = icc;
            if(NSPIN==4) ic = icc/2;
            string orbital_type;
            if (ic == 0)  orbital_type = "S";
            else if (ic == 1) orbital_type = "P";
            else if (ic == 2) orbital_type = "D";
			else if (ic == 3) orbital_type = "F";//mohan add 2009-12-15
            else
            {
				ofs_warning << "\n nchi = " << ucell.atoms[it].nchi << endl;
                WARNING_QUIT("WF_atomic::print_PAOs", "unknown PAO type.");
            }

            stringstream ss;
            ss << global_out_dir
            << ucell.atoms[it].label << "/"
            << ucell.atoms[it].label << "-"
            << orbital_type << ".ORBITAL";

            ofstream ofs(ss.str().c_str());
            ofs << "Mesh " << ucell.atoms[it].msh;
            ofs << "\n" << setw(15) << "Radial"
            << setw(15) << "Psi"
            << setw(15) << "Rab";

            for (int i=0;i<ucell.atoms[it].msh;i++)
            {
                ofs << "\n" << setw(15) << ucell.atoms[it].r[i]
                << setw(15) << ucell.atoms[it].chi(icc,i)
                << setw(15) << ucell.atoms[it].rab[i];
            }
            ofs.close();
        }
        // end out put
    }
    return;
}


//===================================================================
// This routine computes an estimate of the start_ WF_atomictions
// from superposition of atomic WF_atomictions or random wave functions.
//===================================================================
// from wfcinit.f90

void WF_atomic::check_psi(const ComplexMatrix *evc)const
{
    cout<<"\n Check psi : \n";

    for (int iw=0;iw<NBANDS;iw++)
    {
        double sum_evc = abs2_row(this->evc[0],iw);
        OUT("iw",iw);
        OUT("sum_evc",sum_evc);
    }

    for (int ik=0;ik<kv.nks;ik++)
    {
        double sum_evc = abs2(this->evc[ik]);
        OUT("ik",ik);
        OUT("sum_evc2",sum_evc);
    }
    //QUIT();
}

void WF_atomic::atomic_wfc
(	const int ik,
  const int np,
  const int lmax_wfc,
  ComplexMatrix &wfcatom,
  const realArray &table_q,
  const int &table_dimension,
  const double &dq
)const
{
    if (test_wf>3) TITLE("WF_atomic","atomic_wfc");
    timer::tick("WF_atomic","atomic_wfc");
    //=========================================================
    // This routine computes the superposition of atomic
    // WF_atomictions for a given k-point.
    //=========================================================
    const int total_lm = (lmax_wfc + 1) * (lmax_wfc + 1);
    matrix ylm(total_lm, np);
    complex<double> *aux = new complex<double>[np];
    double *chiaux = new double[1];

    Vector3<double> *gk = new Vector3 <double> [np];
    for (int ig=0;ig<np;ig++)
    {
        gk[ig] = WF_atomic::get_1qvec_cartesian(ik, ig);
    }
    //ylm = spherical harmonics functions
    Mathzone::Ylm_Real(total_lm, np, gk, ylm);
    int index = 0;
    //---------------------------------------------------------
    // Calculate G space 3D wave functions
    //---------------------------------------------------------
    double *flq = new double[np];
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
			//cout << "\n it = " << it << " ia = " << ia << endl;
            complex<double> *sk = this->get_sk(ik, it, ia);
            //-------------------------------------------------------
            // Calculate G space 3D wave functions
            //-------------------------------------------------------
            for (int iw = 0;iw < ucell.atoms[it].nchi;iw++)
            {
                if (ucell.atoms[it].oc[iw] >= 0.0)
                {
                    const int l = ucell.atoms[it].lchi[iw];
                    complex<double> lphase = pow(NEG_IMAG_UNIT, l);
                    //-----------------------------------------------------
                    //  the factor i^l MUST BE PRESENT in order to produce
                    //  WF_atomictions for k=0 that are real in real space
                    //-----------------------------------------------------

                      //---------------------------------------------------------
                      // flq = radial fourier transform of atomic orbitals chi
                      //---------------------------------------------------------
                    for (int ig=0; ig<np; ig++)
                    {
                        flq[ig] =
                            Mathzone::Polynomial_Interpolation(table_q,
                                                               it, iw, table_dimension, dq, gk[ig].norm() * ucell.tpiba );
                    }

                    if(NSPIN==4)
                    {
                        if(ucell.atoms[it].has_so)
                        {
                            const double j = ucell.atoms[it].jchi[iw];
                            if ( !(DOMAG||DOMAG_Z))
                            {//atomic_wfc_so
                                double fact[2];
                                for(int m=-l-1;m<l+1;m++)
                                {
                                   fact[0] = soc.spinor(l,j,m,0);
                                   fact[1] = soc.spinor(l,j,m,1);
                                   if (fabs(fact[0])>1e-8||fabs(fact[1])>1e-8)
                                   {
                                      for(int is=0;is<2;is++)
                                      {
                                          if(fabs(fact[is])>1e-8)
                                          {
                                              const int ind = ppcell.lmaxkb + soc.sph_ind(l,j,m,is);
                                              ZEROS(aux, np);
                                              for(int n1=0;n1<2*l+1;n1++){
                                                 const int lm = l*l +n1;
                                                 if(fabs(soc.rotylm(n1,ind))>1e-8)
                                                   for(int ig=0; ig<np;ig++) 
                                                      aux[ig] += soc.rotylm(n1,ind)* ylm(lm,ig);
                                              }
                                              for(int ig=0; ig<np;ig++)
                                                 wfcatom(index, ig + wf.npwx*is ) = lphase * fact[is] * sk[ig] * aux[ig] * flq[ig];
                                          }
                                          else 
                                            for(int ig=0; ig<np;ig++) wfcatom(index,ig+ wf.npwx*is) = complex<double>(0.0 , 0.0);
                                      }//is
                                      index++;
                                   }//if
                                }//m
                            }//if
                            else
                            {//atomic_wfc_so_mag

                              double alpha, gamma;
                              complex<double> fup,fdown;
                              int nc;
                              //This routine creates two functions only in the case j=l+1/2 or exit in the other case  
                              if(fabs(j-l+0.5<1e-4)) continue;
                              delete[] chiaux;
                              chiaux = new double [np];
                              //Find the functions j= l- 1/2
                              if(l==0) 
                                 for(int ig=0;ig<np;ig++){
                                    chiaux[ig] = flq[ig];
                                 }
                              else
                              {
                                 for(int ib = 0;ib < ucell.atoms[it].nchi;ib++)
                                 {
                                    if((ucell.atoms[it].lchi[ib] == l)&&(fabs(ucell.atoms[it].jchi[ib]-l+0.5)<1e-4))
                                    {
                                       nc=ib;
                                       break;
                                    }
                                 }
                                 for(int ig=0;ig<np;ig++)
                                 {//Average the two functions
                                    chiaux[ig] =  l * 
                                         Mathzone::Polynomial_Interpolation(table_q,
                                                               it, nc, table_dimension, dq, gk[ig].norm() * ucell.tpiba );

                                    chiaux[ig] += flq[ig] * (l+1.0) ;
                                    chiaux[ig] *= 1/(2.0*l+1.0);
                                 }
                              }
                              //and construct the starting wavefunctions as in the noncollinear case.
                              alpha = soc.angle1[it];
                              gamma = -1 * soc.angle2[it] + 0.5 * PI;

                              for(int m = 0;m<2*l+1;m++)
                              {
                                 const int lm = l*l +m;
                                 if(index+2*l+1>ucell.natomwfc) WARNING_QUIT("wf.atomic_wfc()","error: too many wfcs");
                                 for(int ig = 0;ig<np;ig++)
                                 {
                                     aux[ig] = sk[ig] * ylm(lm,ig) * chiaux[ig];
                                 }
                                 //rotate wfc as needed
                                 //first rotation with angle alpha around (OX)
                                 for(int ig = 0;ig<np;ig++)
                                 {
                                     fup = cos(0.5 * alpha) * aux[ig];
                                     fdown = IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                     //build the orthogonal wfc
                                     //first rotation with angle (alpha + PI) around (OX)
                                     wfcatom(index,ig) = (cos(0.5 * gamma) + IMAG_UNIT * sin(0.5*gamma)) * fup;
                                     wfcatom(index,ig+ wf.npwx) = (cos(0.5 * gamma) - IMAG_UNIT * sin(0.5*gamma)) * fdown;
                                     //second rotation with angle gamma around(OZ)
                                     fup = cos(0.5 * (alpha + PI))*aux[ig];
                                     fdown = IMAG_UNIT * sin(0.5 * (alpha + PI))*aux[ig];
                                     wfcatom(index+2*l+1,ig) = (cos(0.5*gamma) + IMAG_UNIT*sin(0.5*gamma))*fup;
                                     wfcatom(index+2*l+1,ig+ wf.npwx) = (cos(0.5*gamma) - IMAG_UNIT*sin(0.5*gamma))*fdown;
                                 }
                                 index++;
                              }
                              index += 2*l +1;
                            }
                        }
                        else 
                        {//atomic_wfc_nc
                            double alpha, gamman;
                            complex<double> fup, fdown;
                            alpha = soc.angle1[it];
                            gamman = -soc.angle2[it] + 0.5*PI;
                            for(int m = 0;m<2*l+1;m++)
                            {
                                const int lm = l*l +m;
                                if(index+2*l+1>ucell.natomwfc) WARNING_QUIT("wf.atomic_wfc()","error: too many wfcs");
                                for(int ig = 0;ig<np;ig++)
                                {
                                     aux[ig] = sk[ig] * ylm(lm,ig) * flq[ig];
                                }
                                //rotate function
                                //first, rotation with angle alpha around(OX)
                                for(int ig = 0;ig<np;ig++)
                                {
                                     fup = cos(0.5*alpha) * aux[ig];
                                     fdown = IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                     //build the orthogonal wfc
                                     //first rotation with angle(alpha+PI) around(OX)
                                     wfcatom(index,ig) = (cos(0.5 * gamman) + IMAG_UNIT * sin(0.5*gamman)) * fup;
                                     wfcatom(index,ig+ wf.npwx) = (cos(0.5 * gamman) - IMAG_UNIT * sin(0.5*gamman)) * fdown;
                                     //second rotation with angle gamma around(OZ)
                                     fup = cos(0.5 * (alpha + PI)) * aux[ig];
                                     fdown = IMAG_UNIT * sin(0.5 * (alpha + PI)) * aux[ig];
                                     wfcatom(index+2*l+1,ig) = (cos(0.5*gamman) + IMAG_UNIT*sin(0.5*gamman))*fup;
                                     wfcatom(index+2*l+1,ig+ wf.npwx) = (cos(0.5*gamman) - IMAG_UNIT*sin(0.5*gamman))*fdown;
                                }
                                index++;
                            }
                            index += 2*l+1;
                        }
                       // index++;
                    }
                    else{//LSDA and nomagnet case 

                      for (int m=0; m<2*l+1; m++)
                      {
                        const int lm = l * l + m;
                        for (int ig=0;ig<np;ig++)
                        {
                            wfcatom(index, ig) = lphase * sk [ig] * ylm(lm, ig) * flq[ig];
                            // very useful for debug, don't delete this
//							if(i==24 && index==0){
//							cout << "\n wfcatom(" << index<<","<<i<<") = " << wfcatom(index,i);
//							cout << "\n ylm(" << lm <<","<<i<<") = " << ylm (lm,i)
//							    << " sk[i]=" <<sk[i] << " chiq=" << chiq (it,iw,i) ;
//							}
                        }
                        index++;
                      }//end m
                    }
                } //end if
            } //end iw
			delete [] sk;
        } //end ia //mohan modify 2007-11-7
    } // end nt

	if(test_wf)OUT(ofs_running,"wf_index",index);

    if (index != ucell.natomwfc)
    {
        WARNING_QUIT("wf.atomic_wfc()","index != ucell.natomwfc");
    }
    delete[] flq;
    delete [] gk;
    delete [] aux;
    delete[] chiaux;
    timer::tick("WF_atomic","atomic_wfc");
    return;
} //end subroutine atomic_wfc

void WF_atomic::random(ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const
{
    assert(iw_start >= 0);
    assert(psi.nr >= iw_end);
    const int ng = kv.ngk[ik];
    for (int iw = iw_start ;iw < iw_end;iw++)
    {
        for (int ig = 0;ig < ng;ig++)
        {
            const double rr = std::rand();
            const double arg= TWO_PI * std::rand();
            Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig)];
            psi(iw,ig) = complex<double>(rr * cos(arg), rr * sin(arg)) / (v3 * v3 + 1.0);
        }
        if(NPOL==2)for (int ig = wf.npwx;ig < wf.npwx + ng;ig++)
        {
            const double rr = std::rand();
            const double arg= TWO_PI * std::rand();
            Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig-wf.npwx)];
            psi(iw,ig) = complex<double>(rr * cos(arg), rr * sin(arg)) / (v3 * v3 + 1.0);
        }
    }
    return;
}

