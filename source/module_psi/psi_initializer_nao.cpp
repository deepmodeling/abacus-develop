#include "psi_initializer_nao.h"
#include <fstream>
// numerical algorithm support
#include "module_base/math_integral.h" // for numerical integration
// numerical algorithm support
#include "module_base/math_polyint.h" // for polynomial interpolation
#include "module_base/math_ylmreal.h" // for real spherical harmonics
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
// three global variables definition
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

psi_initializer_nao::psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in) : psi_initializer(sf_in, pw_wfc_in)
{
	this->set_method("nao");
    // find correct dimension for ovlp_flzjlq
    int dim1 = GlobalC::ucell.ntype;
    int dim2 = 0; // dim2 should be the maximum number of zeta for each atomtype
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        int nzeta = 0;
        for (int l = 0; l < GlobalC::ucell.atoms[it].nwl+1; l++)
        {
            nzeta += GlobalC::ucell.atoms[it].l_nchi[l];
        }
        dim2 = (nzeta > dim2) ? nzeta : dim2;
    }
    if (dim2 == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_nao::psi_initializer_nao", "there is not ANY numerical atomic orbital read in present system, quit.");
    }
    int dim3 = GlobalV::NQX;
    // allocate memory for ovlp_flzjlq
    this->ovlp_flzjlq.create(dim1, dim2, dim3);
    this->ovlp_flzjlq.zero_out();
}

psi_initializer_nao::~psi_initializer_nao() {}


void psi_initializer_nao::set_orbital_files(std::string* orbital_files)
{
	ModuleBase::timer::tick("psi_initializer_nao", "set_orbital_files");
	for (int itype = 0; itype < GlobalC::ucell.ntype; itype++)
    {
		this->orbital_files.push_back(orbital_files[itype]);
	}
	ModuleBase::timer::tick("psi_initializer_nao", "set_orbital_files");
}


void psi_initializer_nao::cal_ovlp_flzjlq()
{
	ModuleBase::timer::tick("psi_initializer_nao", "cal_ovlp_flzjlq");
    this->ovlp_flzjlq.zero_out();
	for(int it=0; it<GlobalC::ucell.ntype; it++)
	{
		int ic=0;
		for(int l=0; l<GlobalC::ucell.atoms[it].nwl+1; l++)
		{
			for(int izeta=0; izeta<GlobalC::ucell.atoms[it].l_nchi[l]; izeta++)
			{
				double* ovlp_flzjlq_q = new double[GlobalV::NQX];
                double* qgrid = new double[GlobalV::NQX];
                for (int iq = 0; iq < GlobalV::NQX; iq++)
                {
                    qgrid[iq] = iq*GlobalV::DQ;
                }
                this->sbt.direct(l, 
								 GlobalC::ucell.atoms[it].n_rgrid[ic], 
								 GlobalC::ucell.atoms[it].rgrid[ic], 
								 GlobalC::ucell.atoms[it].flz[ic], 
								 GlobalV::NQX, qgrid, ovlp_flzjlq_q);
				for(int iq = 0; iq < GlobalV::NQX; iq++)
				{
					this->ovlp_flzjlq(it, ic, iq) = ovlp_flzjlq_q[iq];
				}
				delete[] ovlp_flzjlq_q;
				++ic;
			}
		}
	}
	if(GlobalV::MY_RANK==0)
	{
		for(int it = 0; it < GlobalC::ucell.ntype; it++)
		{
			std::stringstream ss;
			ss<<GlobalV::global_out_dir<<GlobalC::ucell.atoms[it].label<< "/LOCAL_G.dat";
			std::ofstream ofs(ss.str().c_str());
			for(int iq = 0; iq < GlobalV::NQX; iq++)
			{
				int ic=0;
				double energy_q = pow((double)iq*GlobalV::DQ, 2);
				ofs<<energy_q*GlobalC::ucell.tpiba2;
				for(int l = 0; l<GlobalC::ucell.atoms[it].nwl + 1; l++)
				{
					for(int N=0; N<GlobalC::ucell.atoms[it].l_nchi[l]; N++)
					{
						ofs<<" "<<ovlp_flzjlq(it, ic, iq);
						++ic;
					}
				}
				ofs<<std::endl;
			}
			ofs.close();
		}
	}
	ModuleBase::timer::tick("psi_initializer_nao", "cal_ovlp_flzjlq");
}

psi::Psi<std::complex<double>>* psi_initializer_nao::cal_psig(int ik)
{
	ModuleBase::timer::tick("psi_initializer_nao", "initialize");
	assert(ik>=0);
	this->psig->fix_k(ik);
	const int npw = this->pw_wfc->npwk[ik];
	const int total_lm = ( GlobalC::ucell.lmax + 1) * ( GlobalC::ucell.lmax + 1);
	ModuleBase::matrix ylm(total_lm, npw);
	std::complex<double> *aux = new std::complex<double>[npw];
	double *chiaux = nullptr;

	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for(int ig=0;ig<npw;ig++)
	{
		gk[ig] = this->pw_wfc->getgpluskcar(ik, ig);
	}

	ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);

	//int index = 0;
	double *ovlp_flzjlg = new double[npw];
	int ibasis=0;
	for (int it = 0; it < GlobalC::ucell.ntype; it++)
	{
		for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
		{
/* HERE LOOP OVER ALL ATOMIS */
            std::complex<double>* sk = this->sf->get_sk(ik, it, ia, this->pw_wfc);
            int ic = 0; // ic is a flatten index of chi, therefore it is defined here.
            for(int L = 0; L < GlobalC::ucell.atoms[it].nwl+1; L++)
			{
				std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, L); //mohan 2010-04-19
				for(int N=0; N < GlobalC::ucell.atoms[it].l_nchi[L]; N++)
				{
/* HERE LOOP OVER ALL NAOS */
					for(int ig=0; ig<npw; ig++)
					{
						ovlp_flzjlg[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(this->ovlp_flzjlq,
						it, ic, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * GlobalC::ucell.tpiba );
					}
/* FOR EVERY NAO IN EACH ATOM */
					if(GlobalV::NSPIN==4)
					{
/* FOR EACH SPIN CHANNEL */
						for(int is_N = 0; is_N < 2; is_N++)
						//for(int is_N = 0; is_N < 1; is_N++)
						{
							if(L==0 && is_N==1) continue;
							if(GlobalC::ucell.atoms[it].ncpp.has_so)
							{
								const double j = std::abs(double(L+is_N) - 0.5);
								if (!(GlobalV::DOMAG||GlobalV::DOMAG_Z))
								{//atomic_wfc_so
									for(int m=0; m<2*L+1; m++)
									{
										std::cout<<"ibasis: "<<ibasis<<std::endl;
										const int lm = L*L+m;
										for(int ig=0; ig<npw; ig++)
										{
											//if(is_N==0)
											(*(this->psig))(ibasis, ig) =
											lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig];
											//else
                                            (*(this->psig))(ibasis + 1, ig + this->pw_wfc->npwk_max)
                                                = lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig];
                                        }
                                        ibasis += 2;
                                    }
								}//if
								else
								{//atomic_wfc_so_mag
									double alpha, gamma;
									std::complex<double> fup,fdown;
                              		//int nc;
                              		//This routine creates two functions only in the case j=l+1/2 or exit in the other case
									if(fabs(j - L + 0.5) < 1e-4) continue;
									delete[] chiaux;
									chiaux = new double [npw];
                              		//Find the functions j= l- 1/2
									if(L==0)
									for(int ig=0;ig<npw;ig++){
										chiaux[ig] = ovlp_flzjlg[ig];
									}
									else
									{
										for(int ig=0;ig<npw;ig++)
										{//Average the two functions
											chiaux[ig] =  L *
												ModuleBase::PolyInt::Polynomial_Interpolation(this->ovlp_flzjlq,
												it, ic, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * GlobalC::ucell.tpiba );

											chiaux[ig] += ovlp_flzjlg[ig] * (L+1.0) ;
											chiaux[ig] *= 1/(2.0*L+1.0);
										}
									}
									alpha = GlobalC::ucell.atoms[it].angle1[ia];
									gamma = -1 * GlobalC::ucell.atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
									for(int m = 0;m<2*L+1;m++)
									{
										const int lm = L*L +m;
                                        for (int ig = 0; ig < npw; ig++)
                                        {
                                            aux[ig] = sk[ig] * ylm(lm,ig) * chiaux[ig];
                                        }
										for(int ig = 0;ig<npw;ig++)
										{
											fup = cos(0.5 * alpha) * aux[ig];
											fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
											//build the orthogonal wfc
											//first rotation with angle (alpha + ModuleBase::PI) around (OX)
											(*(this->psig))(ibasis,ig) = (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5*gamma)) * fup;
                                            (*(this->psig))(ibasis, ig + this->pw_wfc->npwk_max)
                                                = (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown;
                                            // second rotation with angle gamma around(OZ)
                                            fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                            fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI))*aux[ig];
											(*(this->psig))(ibasis+2*L+1,ig) = (cos(0.5*gamma) + ModuleBase::IMAG_UNIT*sin(0.5*gamma))*fup;
                                            (*(this->psig))(ibasis + 2 * L + 1, ig + this->pw_wfc->npwk_max)
                                                = (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown;
                                        }
                                        ibasis++;
                                    }
									ibasis += 2*L +1;
								} // end else INPUT.starting_spin_angle || !GlobalV::DOMAG
							} // end if GlobalC::ucell.atoms[it].has_so
							else
							{//atomic_wfc_nc
								double alpha, gamman;
								std::complex<double> fup, fdown;
								//alpha = GlobalC::ucell.magnet.angle1_[it];
								//gamman = -GlobalC::ucell.magnet.angle2_[it] + 0.5*ModuleBase::PI;
								alpha = GlobalC::ucell.atoms[it].angle1[ia];
								gamman = -GlobalC::ucell.atoms[it].angle2[ia] + 0.5*ModuleBase::PI;
								for(int m = 0; m < 2*L+1; m++)
								{
									const int lm = L*L +m;
                                    for (int ig = 0; ig < npw; ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm,ig) * ovlp_flzjlg[ig];
                                    }
                                    //rotate function
									//first, rotation with angle alpha around(OX)
									for(int ig = 0; ig < npw; ig++)
									{
										fup = cos(0.5*alpha) * aux[ig];
										fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
										//build the orthogonal wfc
										//first rotation with angle(alpha+ModuleBase::PI) around(OX)
										(*(this->psig))(ibasis,ig) = (cos(0.5 * gamman) + ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fup;
                                        (*(this->psig))(ibasis, ig + this->pw_wfc->npwk_max)
                                            = (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown;
                                        // second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
										(*(this->psig))(ibasis+2*L+1,ig) = (cos(0.5*gamman) + ModuleBase::IMAG_UNIT*sin(0.5*gamman))*fup;
                                        (*(this->psig))(ibasis + 2 * L + 1, ig + this->pw_wfc->npwk_max)
                                            = (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown;
                                    } // end ig
                                    ibasis++;
                                } // end m
								ibasis += 2*L+1;
							} // end else GlobalC::ucell.atoms[it].has_so
						} // end for is_N
                    } // end if GlobalV::NONCOLIN
                    else{//LSDA and nomagnet case
	/* DOES NOT DISTINGUISH m QUANTUM NUMBER FOR CHI */
						for(int m = 0; m < 2*L+1; m++)
						{
							const int lm = L*L+m;
							for(int ig=0; ig<npw; ig++)
							{
								(*(this->psig))(ibasis, ig) =
								lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig];
							}
							++ibasis;
						}
					}
					++ic;
				} // end for N
			} // end for L
			delete[] sk;
		} // end for ia
	} // end for it
	delete[] ovlp_flzjlg;
	delete[] aux;
	delete[] chiaux;
	delete[] gk;
	/* complement the rest of bands if there are */
	if(this->get_nbands_complem() > 0)
	{
		this->random_t(this->psig->get_pointer(), ibasis, this->psig->get_nbands(), ik, this->pw_wfc);
	}
	ModuleBase::timer::tick("psi_initializer_nao", "initialize");
	
#ifdef PSI_INITIALIZER_TEST
	this->write_psig();
#endif
	return this->psig;
}