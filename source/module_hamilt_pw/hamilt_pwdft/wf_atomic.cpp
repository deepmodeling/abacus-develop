#include "wf_atomic.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
#include <complex>
#include "module_base/timer.h"
#include "module_base/tool_quit.h"

WF_atomic::WF_atomic()
{
    pw_seed = 0;
}

WF_atomic::~WF_atomic()
{
	if(GlobalV::test_deconstructor)
	{
		std::cout << " ~WF_atomic()" << std::endl;
	}
    if(this->evc!=nullptr)
    {
        delete[] evc;
    }
    if(this->wanf2!= nullptr)
    {
        delete[] wanf2;
    }
    if(this->psi != nullptr)
    {
        delete psi;
    }
}

//==========================================================
// MEMBER FUNCTION :
// NAME : init_at_1(init a table with the radial Fourier
// transform of the atomic WF_atomictions)
//==========================================================
void WF_atomic::init_at_1(void)
{
    if (GlobalV::test_wf) ModuleBase::TITLE("WF_atomic","init_at_1");
    ModuleBase::timer::tick("WF_atomic","init_at_1");

	GlobalV::ofs_running << "\n Make real space PAO into reciprocal space." << std::endl;

    this->print_PAOs();

//----------------------------------------------------------
// EXPLAIN : Find the type of atom that has most mesh points.
//----------------------------------------------------------
    int ndm = 0;
    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {
        ndm = (GlobalC::ucell.atoms[it].ncpp.msh > ndm) ? GlobalC::ucell.atoms[it].ncpp.msh : ndm;
    }
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max mesh points in Pseudopotential",ndm);

    // needed to normalize atomic wfcs (not a bad idea in general and
    // necessary to compute correctly lda+U projections)
    GlobalC::ppcell.tab_at.zero_out();
//----------------------------------------------------------
// EXPLAIN : If use gauss orbitals to represent aotmic
// orbitals (controlled by parameters)
//
// USE GLOBAL CLASS VARIABLES :
// NAME : GlobalC::ucell.atoms.nchi
// NAME : GlobalC::ucell.atoms.msh(number of mesh points)
// NAME : GlobalC::ucell.atoms.r
//----------------------------------------------------------
    const int startq = 0;
    const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);
    double *aux = new double[ndm];
    double *vchi = new double[ndm];

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dq(describe PAO in reciprocal space)",GlobalV::DQ);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max q",GlobalV::NQX);

    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {
		Atom* atom = &GlobalC::ucell.atoms[it];

		GlobalV::ofs_running <<"\n number of pseudo atomic orbitals for "
		<< atom->label << " is " << atom->ncpp.nchi << std::endl;

        for (int ic=0; ic<atom->ncpp.nchi ;ic++)
        {
			//std::cout << "\n T=" << it << " ic=" << ic << std::endl;
            int nmesh;
            if(GlobalV::PSEUDO_MESH)
                nmesh = atom->ncpp.mesh;
            else
                nmesh = atom->ncpp.msh;

            // check the unit condition
            double *inner_part = new double[nmesh];
            for (int ir=0; ir<nmesh; ir++)
            {
                inner_part[ir] = atom->ncpp.chi(ic,ir) * atom->ncpp.chi(ic,ir);
            }
            double unit = 0.0;
            ModuleBase::Integral::Simpson_Integral(nmesh, inner_part, atom->ncpp.rab, unit);
            delete[] inner_part;

			GlobalV::ofs_running << " the unit of pseudo atomic orbital is " << unit;

            //=================================
            // normalize radial wave functions
            //=================================
            for (int ir=0; ir<nmesh; ir++)
            {
                atom->ncpp.chi(ic,ir) /= sqrt(unit);
            }

            //===========
            // recheck
            //===========
            inner_part = new double[nmesh];
            for (int ir=0; ir<nmesh; ir++)
            {
                inner_part[ir] = atom->ncpp.chi(ic,ir) * atom->ncpp.chi(ic,ir);
            }
            unit = 0.0;
            ModuleBase::Integral::Simpson_Integral(nmesh, inner_part, atom->ncpp.rab, unit);
            delete[] inner_part;

			GlobalV::ofs_running << ", renormalize to " << unit << std::endl;

            if (atom->ncpp.oc[ic] >= 0.0)
            {
                const int l = atom->ncpp.lchi[ic];
                for (int iq=startq; iq<GlobalV::NQX; iq++)
                {
                    const double q = GlobalV::DQ * iq;
                    ModuleBase::Sphbes::Spherical_Bessel(atom->ncpp.msh, atom->ncpp.r, q, l, aux);
                    for (int ir = 0;ir < atom->ncpp.msh;ir++)
                    {
                        vchi[ir] = atom->ncpp.chi(ic,ir) * aux[ir] * atom->ncpp.r[ir];
                    }

                    double vqint = 0.0;
                    ModuleBase::Integral::Simpson_Integral(atom->ncpp.msh, vchi, atom->ncpp.rab, vqint);

                    GlobalC::ppcell.tab_at(it, ic, iq) =  vqint * pref;
                    //				if( it == 0 && ic == 0 )
                    //				{
                    //
                    //					for (ir = 0;ir < GlobalC::ucell.atoms[it].msh;ir++)
                    //						GlobalV::ofs_running << std::setprecision(20) << "\n vchi(" << ir << ")=" << vchi[ir];
                    //					GlobalV::ofs_running<<"\n aux[0] = "<<aux[0];
                    //					GlobalV::ofs_running<<"\n msh = "<< GlobalC::ucell.atoms[it].msh;
                    //					GlobalV::ofs_running<<"\n tab_at : "<<ppcell.tab_at(it, ic, iq) ;
                    //					GlobalV::ofs_running<<"\n pref = "<<pref;
                    //				}
                } // enddo
            } // endif
        } // enddo
    } // enddo

    delete [] aux;
    delete [] vchi;
    ModuleBase::timer::tick("WF_atomic","init_at_1");
    return;
}// end init_at_1

void WF_atomic::print_PAOs(void)const
{
    if (GlobalV::MY_RANK!=0) return;
    for (int it=0; it<GlobalC::ucell.ntype; it++)
    {
        for (int icc=0; icc<GlobalC::ucell.atoms[it].ncpp.nchi ;icc++)
        {
            int ic = icc;
            if(GlobalV::NSPIN==4) ic = icc/2;
            std::string orbital_type;
            if (ic == 0)  orbital_type = "S";
            else if (ic == 1) orbital_type = "P";
            else if (ic == 2) orbital_type = "D";
			else if (ic == 3) orbital_type = "F"; // mohan add 2009-12-15
			else if (ic == 4) orbital_type = "G"; // liuyu add 2021-05-07
            else
            {
				GlobalV::ofs_warning << "\n nchi = " << GlobalC::ucell.atoms[it].ncpp.nchi << std::endl;
                ModuleBase::WARNING_QUIT("WF_atomic::print_PAOs", "unknown PAO type.");
            }

            std::stringstream ss;
            ss << GlobalV::global_out_dir
            << GlobalC::ucell.atoms[it].label << "/"
            << GlobalC::ucell.atoms[it].label << "-"
            << orbital_type << ".ORBITAL";

            std::ofstream ofs(ss.str().c_str());
            ofs << "Mesh " << GlobalC::ucell.atoms[it].ncpp.msh;
            ofs << "\n" << std::setw(15) << "Radial"
            << std::setw(15) << "Psi"
            << std::setw(15) << "Rab";

            for (int i=0;i<GlobalC::ucell.atoms[it].ncpp.msh;i++)
            {
                ofs << "\n" << std::setw(15) << GlobalC::ucell.atoms[it].ncpp.r[i]
                << std::setw(15) << GlobalC::ucell.atoms[it].ncpp.chi(icc,i)
                << std::setw(15) << GlobalC::ucell.atoms[it].ncpp.rab[i];
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

void WF_atomic::check_evc()const
{
    std::cout<<"\n Check psi : \n";

    for (int iw=0;iw<GlobalV::NBANDS;iw++)
    {
        double sum_evc = abs2_row(this->evc[0],iw);
        ModuleBase::GlobalFunc::OUT("iw",iw);
        ModuleBase::GlobalFunc::OUT("sum_evc",sum_evc);
    }

    for (int ik=0;ik<GlobalC::kv.nks;ik++)
    {
        double sum_evc = abs2(this->evc[ik]);
        ModuleBase::GlobalFunc::OUT("ik",ik);
        ModuleBase::GlobalFunc::OUT("sum_evc2",sum_evc);
    }
    //ModuleBase::QUIT();
}

void WF_atomic::atomic_wfc
(	const int ik,
  const int np,
  const int lmax_wfc,
  ModuleBase::ComplexMatrix &wfcatom,
  const ModuleBase::realArray &table_q,
  const int &table_dimension,
  const double &dq
)const
{
    if (GlobalV::test_wf>3) ModuleBase::TITLE("WF_atomic","atomic_wfc");
    ModuleBase::timer::tick("WF_atomic","atomic_wfc");
    //=========================================================
    // This routine computes the superposition of atomic
    // WF_atomictions for a given k-point.
    //=========================================================
    const int total_lm = (lmax_wfc + 1) * (lmax_wfc + 1);
    ModuleBase::matrix ylm(total_lm, np);
    std::complex<double> *aux = new std::complex<double>[np];
    double *chiaux = nullptr;

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3 <double> [np];
    for (int ig=0;ig<np;ig++)
    {
        gk[ig] = WF_atomic::get_1qvec_cartesian(ik, ig);
    }
    //ylm = spherical harmonics functions
    ModuleBase::YlmReal::Ylm_Real(total_lm, np, gk, ylm);
    int index = 0;
    //---------------------------------------------------------
    // Calculate G space 3D wave functions
    //---------------------------------------------------------
    double *flq = new double[np];
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
			//std::cout << "\n it = " << it << " ia = " << ia << std::endl;
            std::complex<double> *sk = this->get_sk(ik, it, ia,GlobalC::wfcpw);
            //-------------------------------------------------------
            // Calculate G space 3D wave functions
            //-------------------------------------------------------
            for (int iw = 0;iw < GlobalC::ucell.atoms[it].ncpp.nchi;iw++)
            {
                if (GlobalC::ucell.atoms[it].ncpp.oc[iw] >= 0.0)
                {
                    const int l = GlobalC::ucell.atoms[it].ncpp.lchi[iw];
                    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, l);
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
                            ModuleBase::PolyInt::Polynomial_Interpolation(table_q,
                                                               it, iw, table_dimension, dq, gk[ig].norm() * GlobalC::ucell.tpiba );
                    }

                    if(GlobalV::NSPIN==4)
                    {
                        if(GlobalC::ucell.atoms[it].ncpp.has_so)
                        {
                            Soc soc;
						    soc.rot_ylm(l+1);
                            const double j = GlobalC::ucell.atoms[it].ncpp.jchi[iw];
                            if ( !(GlobalV::DOMAG||GlobalV::DOMAG_Z))
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
                                              const int ind = GlobalC::ppcell.lmaxkb + soc.sph_ind(l,j,m,is);
                                              ModuleBase::GlobalFunc::ZEROS(aux, np);
                                              for(int n1=0;n1<2*l+1;n1++){
                                                 const int lm = l*l +n1;
                                                 if(abs(soc.rotylm(n1,ind))>1e-8)
                                                   for(int ig=0; ig<np;ig++)
                                                      aux[ig] += soc.rotylm(n1,ind)* ylm(lm,ig);
                                              }
                                              for(int ig=0; ig<np;ig++)
                                                 wfcatom(index, ig + this->npwx*is ) = lphase * fact[is] * sk[ig] * aux[ig] * flq[ig];
                                          }
                                          else
                                            for(int ig=0; ig<np;ig++) wfcatom(index,ig+ this->npwx*is) = std::complex<double>(0.0 , 0.0);
                                      }//is
                                      index++;
                                   }//if
                                }//m
                            }//if
                            else
                            {//atomic_wfc_so_mag

                                double alpha, gamma;
                                std::complex<double> fup,fdown;
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
                                    for(int ib = 0;ib < GlobalC::ucell.atoms[it].ncpp.nchi;ib++)
                                    {
                                        if((GlobalC::ucell.atoms[it].ncpp.lchi[ib] == l)&&(fabs(GlobalC::ucell.atoms[it].ncpp.jchi[ib]-l+0.5)<1e-4))
                                        {
                                        nc=ib;
                                        break;
                                        }
                                    }
                                    for(int ig=0;ig<np;ig++)
                                    {//Average the two functions
                                        chiaux[ig] =  l *
                                            ModuleBase::PolyInt::Polynomial_Interpolation(table_q,
                                                                it, nc, table_dimension, dq, gk[ig].norm() * GlobalC::ucell.tpiba );

                                        chiaux[ig] += flq[ig] * (l+1.0) ;
                                        chiaux[ig] *= 1/(2.0*l+1.0);
                                    }
                                }
                                //and construct the starting wavefunctions as in the noncollinear case.
                                //alpha = GlobalC::ucell.magnet.angle1_[it];
                                //gamma = -1 * GlobalC::ucell.magnet.angle2_[it] + 0.5 * ModuleBase::PI;
                                alpha = GlobalC::ucell.atoms[it].angle1[ia];
                                gamma = -1 * GlobalC::ucell.atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;

                                for(int m = 0;m<2*l+1;m++)
                                {
                                    const int lm = l*l +m;
                                    if(index+2*l+1>GlobalC::ucell.natomwfc) ModuleBase::WARNING_QUIT("GlobalC::wf.atomic_wfc()","error: too many wfcs");
                                    for(int ig = 0;ig<np;ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm,ig) * chiaux[ig];
                                    }
                                    //rotate wfc as needed
                                    //first rotation with angle alpha around (OX)
                                    for(int ig = 0;ig<np;ig++)
                                    {
                                        fup = cos(0.5 * alpha) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                        //build the orthogonal wfc
                                        //first rotation with angle (alpha + ModuleBase::PI) around (OX)
                                        wfcatom(index,ig) = (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5*gamma)) * fup;
                                        wfcatom(index,ig+ this->npwx) = (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5*gamma)) * fdown;
                                        //second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI))*aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI))*aux[ig];
                                        wfcatom(index+2*l+1,ig) = (cos(0.5*gamma) + ModuleBase::IMAG_UNIT*sin(0.5*gamma))*fup;
                                        wfcatom(index+2*l+1,ig+ this->npwx) = (cos(0.5*gamma) - ModuleBase::IMAG_UNIT*sin(0.5*gamma))*fdown;
                                    }
                                    index++;
                                }
                                index += 2*l +1;
                            }
                        }
                        else
                        {//atomic_wfc_nc
                            double alpha, gamman;
                            std::complex<double> fup, fdown;
                            //alpha = GlobalC::ucell.magnet.angle1_[it];
                            //gamman = -GlobalC::ucell.magnet.angle2_[it] + 0.5*ModuleBase::PI;
                            alpha = GlobalC::ucell.atoms[it].angle1[ia];
                            gamman = -1 * GlobalC::ucell.atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
                            for(int m = 0;m<2*l+1;m++)
                            {
                                const int lm = l*l +m;
                                if(index+2*l+1>GlobalC::ucell.natomwfc) ModuleBase::WARNING_QUIT("GlobalC::wf.atomic_wfc()","error: too many wfcs");
                                for(int ig = 0;ig<np;ig++)
                                {
                                     aux[ig] = sk[ig] * ylm(lm,ig) * flq[ig];
                                }
                                //rotate function
                                //first, rotation with angle alpha around(OX)
                                for(int ig = 0;ig<np;ig++)
                                {
                                     fup = cos(0.5*alpha) * aux[ig];
                                     fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                     //build the orthogonal wfc
                                     //first rotation with angle(alpha+ModuleBase::PI) around(OX)
                                     wfcatom(index,ig) = (cos(0.5 * gamman) + ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fup;
                                     wfcatom(index,ig+ this->npwx) = (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fdown;
                                     //second rotation with angle gamma around(OZ)
                                     fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                     fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                     wfcatom(index+2*l+1,ig) = (cos(0.5*gamman) + ModuleBase::IMAG_UNIT*sin(0.5*gamman))*fup;
                                     wfcatom(index+2*l+1,ig+ this->npwx) = (cos(0.5*gamman) - ModuleBase::IMAG_UNIT*sin(0.5*gamman))*fdown;
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
//							std::cout << "\n wfcatom(" << index<<","<<i<<") = " << wfcatom(index,i);
//							std::cout << "\n ylm(" << lm <<","<<i<<") = " << ylm (lm,i)
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

	if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"wf_index",index);

    if (index != GlobalC::ucell.natomwfc)
    {
        ModuleBase::WARNING_QUIT("GlobalC::wf.atomic_wfc()","index != GlobalC::ucell.natomwfc");
    }
    delete[] flq;
    delete [] gk;
    delete [] aux;
    delete[] chiaux;
    ModuleBase::timer::tick("WF_atomic","atomic_wfc");
    return;
} //end subroutine atomic_wfc

#ifdef __MPI
void WF_atomic::stick_to_pool(double *stick, const int &ir, double *out, ModulePW::PW_Basis_K* wfc_basis) const
{	
	MPI_Status ierror;
    const int is = this->irindex[ir];
	const int ip = wfc_basis->fftixy2ip[ir];
    const int nz = wfc_basis->nz;

	if(ip == 0 && GlobalV::RANK_IN_POOL ==0)
	{
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(ip == GlobalV::RANK_IN_POOL )
	{
		MPI_Recv(stick, nz, MPI_DOUBLE, 0, ir, POOL_WORLD,&ierror);
		for(int iz=0; iz<nz; iz++)
		{
			out[is*nz+iz] = stick[iz];
		}
	}
	else if(GlobalV::RANK_IN_POOL==0)
	{
		MPI_Send(stick, nz, MPI_DOUBLE, ip, ir, POOL_WORLD);
	}

	return;	
}
#endif

void WF_atomic::random(ModuleBase::ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis)
{
    assert(iw_start >= 0);
    assert(psi.nr >= iw_end);
    const int ng = GlobalC::kv.ngk[ik];
#ifdef __MPI
// #if ((defined __CUDA) || (defined __ROCM))
    // if(pw_seed > 0)//qianrui add 2021-8-13
    // {
    //     srand(unsigned(pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
    // }
// #else
    if(pw_seed > 0)//qianrui add 2021-8-13
    {
        srand(unsigned(pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
        const int nxy = wfc_basis->fftnxy;
        const int nz = wfc_basis->nz;
        const int nstnz = wfc_basis->nst*nz;

        double *stickrr = new double[nz];
        double *stickarg = new double[nz];
        double *tmprr = new double[nstnz];
        double *tmparg = new double[nstnz];
        for (int iw = iw_start ;iw < iw_end;iw++)
        {
            int startig = 0;
            for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
            {
            
	            for(int ir=0; ir < nxy; ir++)
	            {
                    if(wfc_basis->fftixy2ip[ir] < 0) continue;
	            	if(GlobalV::RANK_IN_POOL==0)
	            	{
	            		for(int iz=0; iz<nz; iz++)
	            		{
	            			stickrr[ iz ] = std::rand()/double(RAND_MAX);
                            stickarg[ iz ] = std::rand()/double(RAND_MAX);
	            		}
	            	}
	            	stick_to_pool(stickrr, ir, tmprr, wfc_basis);
                    stick_to_pool(stickarg, ir, tmparg, wfc_basis);
	            }

                for (int ig = 0;ig < ng;ig++)
                {
                    const double rr = tmprr[wfc_basis->getigl2isz(ik,ig)];
                    const double arg= ModuleBase::TWO_PI * tmparg[wfc_basis->getigl2isz(ik,ig)];
                    const double gk2 = GlobalC::wfcpw->getgk2(ik,ig);
                    psi(iw,ig+startig) = std::complex<double>(rr * cos(arg), rr * sin(arg)) / (gk2 + 1.0);
                }
                startig += npwx;
            }
        }
        delete[] stickrr;
        delete[] stickarg;
        delete[] tmprr;
        delete[] tmparg;
    }
    else
    {
// #endif
#else  // !__MPI
        if(pw_seed > 0)//qianrui add 2021-8-13
        {
            srand(unsigned(pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
        }
#endif // __MPI
        for (int iw = iw_start ;iw < iw_end;iw++)
        {
            for (int ig = 0;ig < ng;ig++)
            {
                const double rr = std::rand()/double(RAND_MAX); //qianrui add RAND_MAX
                const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = GlobalC::wfcpw->getgk2(ik,ig);
                psi(iw,ig) = std::complex<double>(rr * cos(arg), rr * sin(arg)) / (gk2 + 1.0);
            }
            if(GlobalV::NPOL==2)for (int ig = this->npwx;ig < this->npwx + ng;ig++)
            {
                const double rr = std::rand()/double(RAND_MAX);
                const double arg= ModuleBase::TWO_PI * std::rand()/double(RAND_MAX);
                const double gk2 = GlobalC::wfcpw->getgk2(ik,ig-this->npwx);
                psi(iw,ig) = std::complex<double>(rr * cos(arg), rr * sin(arg)) / (gk2 + 1.0);
            }
        }
#ifdef __MPI
// #if ((!defined __CUDA) && (!defined __ROCM))
    }
// #endif // ((!defined __CUDA) && (!defined __ROCM))
#endif // __MPI
}

void WF_atomic::atomicrandom(ModuleBase::ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis)const
{
    assert(iw_start >= 0);
    assert(psi.nr >= iw_end);
    const int ng = GlobalC::kv.ngk[ik];
#ifdef __MPI
    if(pw_seed > 0)//qianrui add 2021-8-13
    {
        srand(unsigned(pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
        const int nxy = wfc_basis->fftnxy;
        const int nz = wfc_basis->nz;
        const int nstnz = wfc_basis->nst*nz;

        double *stickrr = new double[nz];
        double *stickarg = new double[nz];
        double *tmprr = new double[nstnz];
        double *tmparg = new double[nstnz];
        for (int iw = iw_start ;iw < iw_end;iw++)
        {
            int startig = 0;
            for(int ipol = 0 ; ipol < GlobalV::NPOL ; ++ipol)
            {
	            for(int ir=0; ir < nxy; ir++)
	            {
                    if(wfc_basis->fftixy2ip[ir] < 0) continue;
	            	if(GlobalV::RANK_IN_POOL==0)
	            	{
	            		for(int iz=0; iz<nz; iz++)
	            		{
	            			stickrr[ iz ] = std::rand()/double(RAND_MAX);
                            stickarg[ iz ] = std::rand()/double(RAND_MAX);
	            		}
	            	}
	            	stick_to_pool(stickrr, ir, tmprr, wfc_basis);
                    stick_to_pool(stickarg, ir, tmparg, wfc_basis);
	            }

                for (int ig = 0;ig < ng;ig++)
                {
                    const double rr = tmprr[wfc_basis->ig2isz[ig]];
                    const double arg= ModuleBase::TWO_PI * tmparg[wfc_basis->ig2isz[ig]];
                    psi(iw,startig+ig) *= (1.0 + 0.05 * std::complex<double>(rr * cos(arg), rr * sin(arg)));
                }
                startig += npwx;
            }
        }
        delete[] stickrr;
        delete[] stickarg;
        delete[] tmprr;
        delete[] tmparg;
    }
    else
    {
#else
        if(pw_seed > 0)//qianrui add 2021-8-13
        {
            srand(unsigned(pw_seed + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL] + ik));
        }
#endif
        double rr, arg;
		for (int iw = iw_start ;iw < iw_end;iw++)
		{
			int startig = 0;
			for(int ip = 0 ; ip < GlobalV::NPOL; ++ip)
			{
				for(int ig = 0 ; ig < npw ; ++ig)
				{
					rr = rand()/double(RAND_MAX);
					arg = ModuleBase::TWO_PI * rand()/double(RAND_MAX);
					psi(iw,startig+ig) *= (1.0 + 0.05 * std::complex<double>(rr * cos(arg), rr * sin(arg)));
				}
				startig += npwx;
			}
		}
#ifdef __MPI
    }
#endif

    return;
}

void WF_atomic::evc_transform_psi()
{
    if(this->evc==nullptr || this->psi != nullptr)
    {
        //ModuleBase::WARNING_QUIT("WF_atomic","no evc or psi is not nullptr, please check!");
        std::cout<<__FILE__<<__LINE__<<" there is no need to transform!"<<std::endl;
        return;
    }
    this->psi = new psi::Psi<std::complex<double>>(GlobalC::kv.nks, this->evc[0].nr, this->evc[0].nc, GlobalC::kv.ngk.data());
    for(int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        this->psi->fix_k(ik);
        std::complex<double> *ppsi = this->psi->get_pointer();
        std::complex<double> *pevc = this->evc[ik].c;
        for(int index = 0; index < this->evc[0].nc * this->evc[0].nr; ++index)
        {
            ppsi[index] =  pevc[index];
        }
    }
    delete[] this->evc;
    this->evc = nullptr;

    return;
}

void WF_atomic::psi_transform_evc()
{
    if(this->psi==nullptr || this->evc != nullptr)
    {
        ModuleBase::WARNING_QUIT("WF_atomic","no psi, please check!");
    }
    this->evc = new ModuleBase::ComplexMatrix [this->psi->get_nk()];
    for(int ik = 0; ik < this->psi->get_nk(); ++ik)
    {
        this->evc[ik].create(this->psi->get_nbands(), this->psi->get_nbasis());
        this->psi->fix_k(ik);
        std::complex<double> *ppsi = this->psi->get_pointer();
        std::complex<double> *pevc = this->evc[ik].c;
        for(int index = 0; index < this->evc[ik].nc * this->evc[ik].nr; ++index)
        {
            pevc[index] =  ppsi[index];
        }
    }
    delete this->psi;
    this->psi = nullptr;

    return;
}
