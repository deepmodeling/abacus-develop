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
// parallel communication
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

template <typename T, typename Device>
#ifdef __MPI
psi_initializer_nao<T, Device>::psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in)
					: psi_initializer<T, Device>(sf_in, pw_wfc_in, p_ucell_in, p_parakpts_in, random_seed_in)
#else
psi_initializer_nao<T, Device>::psi_initializer_nao(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in)
					: psi_initializer<T, Device>(sf_in, pw_wfc_in, p_ucell_in, random_seed_in)
#endif
{
	this->set_method("nao");
}

template <typename T, typename Device>
psi_initializer_nao<T, Device>::~psi_initializer_nao() {}

/*
I don't know why some variables are distributed while others not... for example the orbital_files...
We need not only read and import, but also distribute here
*/
template <typename T, typename Device>
void psi_initializer_nao<T, Device>::set_orbital_files(std::string* orbital_files)
{
	ModuleBase::timer::tick("psi_initializer_nao", "set_orbital_files");
	#ifdef __MPI
	if(GlobalV::MY_RANK == 0)
	{
	#endif
		for (int itype = 0; itype < this->p_ucell->ntype; itype++)
		{
			this->orbital_files.push_back(orbital_files[itype]);
		}
	#ifdef __MPI
	}
	else
	{
		this->orbital_files.resize(this->p_ucell->ntype);
	}
	Parallel_Common::bcast_string(this->orbital_files.data(), this->p_ucell->ntype);
	#endif
	ModuleBase::timer::tick("psi_initializer_nao", "set_orbital_files");
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::create_ovlp_Xjlq()
{
    // find correct dimension for ovlp_flzjlq
    int dim1 = this->p_ucell->ntype;
    int dim2 = 0; // dim2 should be the maximum number of zeta for each atomtype
    for (int it = 0; it < this->p_ucell->ntype; it++)
    {
        int nzeta = 0;
        for (int l = 0; l < this->p_ucell->atoms[it].nwl+1; l++)
        {
            nzeta += this->p_ucell->atoms[it].l_nchi[l];
        }
        dim2 = (nzeta > dim2) ? nzeta : dim2;
    }
    if (dim2 == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_nao<T, Device>::psi_initializer_nao", "there is not ANY numerical atomic orbital read in present system, quit.");
    }
    int dim3 = GlobalV::NQX;
    // allocate memory for ovlp_flzjlq
    this->ovlp_flzjlq.create(dim1, dim2, dim3);
    this->ovlp_flzjlq.zero_out();
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::read_orbital_files()
{
	ModuleBase::timer::tick("psi_initializer_nao", "read_orbital_files");
	#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
	#endif
		for(int it = 0; it < this->p_ucell->ntype; it++)
		{
			// number of chi per atomtype
			int nchi = 0;
			for(int l = 0; l <= this->p_ucell->atoms[it].nwl; l++)
			{
				nchi += this->p_ucell->atoms[it].l_nchi[l];
			}
			
			std::vector<int> n_rgrid_it;
			std::vector<std::vector<double>> rgrid_it;
			std::vector<std::vector<double>> flz_it;

			std::ifstream ifs_it;
			ifs_it.open(GlobalV::global_orbital_dir+this->orbital_files[it]);

			if(!ifs_it)
			{
				GlobalV::ofs_warning<<"psi_initializer_nao<T, Device>::read_orbital_files: cannot open orbital file: "<<this->orbital_files[it]<<std::endl;
				ModuleBase::WARNING_QUIT("psi_initializer_nao<T, Device>::read_orbital_files", "cannot open orbital file.");
			}
			else
			{
				GlobalV::ofs_running<<"psi_initializer_nao<T, Device>::read_orbital_files: reading orbital file: "<<this->orbital_files[it]<<std::endl;
			}
			ifs_it.close();

			int ichi_overall = 0;
			// check nwl and nchi for each rank
			
			for(int l = 0; l <= this->p_ucell->atoms[it].nwl; l++)
			{
				for(int ichi = 0; ichi < this->p_ucell->atoms[it].l_nchi[l]; ichi++)
				{
					int n_rgrid_ichi;
					std::vector<double> rgrid_ichi;
					std::vector<double> flz_ichi;

					GlobalV::ofs_running<<"-------------------------------------- "<<std::endl;
					GlobalV::ofs_running<<" reading orbital of element "<<this->p_ucell->atoms[it].label<<std::endl
										<<" angular momentum l = "<<l<<std::endl
										<<" index of chi = "<<ichi<<std::endl;

					ifs_it.open(GlobalV::global_orbital_dir+this->orbital_files[it]);
					double dr = 0.0;
					char word[80];
					
					while(ifs_it.good())
					{
						ifs_it>>word;
						if(std::strcmp(word, "END")==0) break;
					}
					ModuleBase::CHECK_NAME(ifs_it, "Mesh");
					ifs_it>>n_rgrid_ichi;
					
					if(n_rgrid_ichi%2 == 0) ++n_rgrid_ichi;
					GlobalV::ofs_running<<" number of radial grid = "<<n_rgrid_ichi<<std::endl;

					ModuleBase::CHECK_NAME(ifs_it, "dr");
					ifs_it>>dr;
					GlobalV::ofs_running<<" dr = "<<dr<<std::endl;

					for(int ir = 0; ir < n_rgrid_ichi; ir++)
					{
						rgrid_ichi.push_back(ir*dr);
					}
					GlobalV::ofs_running<<" maximal radial grid point = "<<rgrid_ichi[n_rgrid_ichi-1]<<" Angstrom"<<std::endl;
					
					std::string title1, title2, title3;
					int it_read, l_read, nchi_read;
					bool find = false;
					while(!find)
					{
						if(ifs_it.eof())
						{
							GlobalV::ofs_warning<<" psi_initializer_nao<T, Device>::read_orbital_files: cannot find orbital of element "<<this->p_ucell->atoms[it].label<<std::endl
												<<" angular momentum l = "<<l<<std::endl
												<<" index of chi = "<<ichi<<std::endl;
						}
						ifs_it>>title1>>title2>>title3;
						assert(title1=="Type");
						ifs_it>>it_read>>l_read>>nchi_read;
						if(l_read==l && nchi_read==ichi)
						{
							for(int ir = 0; ir<n_rgrid_ichi; ir++)
							{
								double flz_ichi_ir;
								ifs_it>>flz_ichi_ir;
								flz_ichi.push_back(flz_ichi_ir);
							}
							find = true;
						}
						else
						{
							double discard;
							for(int ir = 0; ir<n_rgrid_ichi; ir++)
							{
								ifs_it>>discard;
							}
						}
					}
					ifs_it.close();
					n_rgrid_it.push_back(n_rgrid_ichi);
					rgrid_it.push_back(rgrid_ichi);
					flz_it.push_back(flz_ichi);
					++ichi_overall;
				}
			}
			this->n_rgrid.push_back(n_rgrid_it);
			this->rgrid.push_back(rgrid_it);
			this->flz.push_back(flz_it);
			GlobalV::ofs_running<<"-------------------------------------- "<<std::endl;
		}
	#ifdef __MPI
	}
	#endif
	// Broadcast n_rgrid, rgrid, flz
	#ifdef __MPI
	if(GlobalV::MY_RANK!=0)
	{
		this->n_rgrid.resize(this->p_ucell->ntype);
	}
	int nchi[this->p_ucell->ntype];
	if(GlobalV::MY_RANK==0)
	{
		for(int it = 0; it < this->p_ucell->ntype; it++)
		{
			nchi[it] = this->n_rgrid[it].size();
		}
	}
	Parallel_Common::bcast_int(nchi, this->p_ucell->ntype);
	if(GlobalV::MY_RANK!=0)
	{
		this->n_rgrid.resize(this->p_ucell->ntype);
		this->rgrid.resize(this->p_ucell->ntype);
		this->flz.resize(this->p_ucell->ntype);
		for(int it = 0; it < this->p_ucell->ntype; it++)
		{
			this->n_rgrid[it].resize(nchi[it]);
			this->rgrid[it].resize(nchi[it]);
			this->flz[it].resize(nchi[it]);
		}
	}
	for(int it = 0; it < this->p_ucell->ntype; it++)
	{
		Parallel_Common::bcast_int(this->n_rgrid[it].data(), nchi[it]);
	}
	if(GlobalV::MY_RANK!=0)
	{
		for(int it = 0; it < this->p_ucell->ntype; it++)
		{
			for(int ichi = 0; ichi < nchi[it]; ichi++)
			{
				this->rgrid[it][ichi].resize(this->n_rgrid[it][ichi]);
				this->flz[it][ichi].resize(this->n_rgrid[it][ichi]);
			}
		}
	}
	for(int it = 0; it < this->p_ucell->ntype; it++)
	{
		for(int ichi = 0; ichi < nchi[it]; ichi++)
		{
			Parallel_Common::bcast_double(this->rgrid[it][ichi].data(), this->n_rgrid[it][ichi]);
			Parallel_Common::bcast_double(this->flz[it][ichi].data(), this->n_rgrid[it][ichi]);
		}
	}
	#endif
	ModuleBase::timer::tick("psi_initializer_nao", "read_orbital_files");
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in)
{
	ModuleBase::timer::tick("psi_initializer_nao", "initialize_only_once");
	this->create_ovlp_Xjlq();
	this->read_orbital_files();
	//this->cal_ovlp_flzjlq(); //because GlobalV::NQX will change during vcrelax, so it should be called in both init and init_after_vc
	ModuleBase::timer::tick("psi_initializer_nao", "initialize_only_once");
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::cal_ovlp_flzjlq()
{
	ModuleBase::timer::tick("psi_initializer_nao", "cal_ovlp_flzjlq");
	//this->read_orbital_files();
    this->ovlp_flzjlq.zero_out();
	for(int it=0; it<this->p_ucell->ntype; it++)
	{
		int ic=0;
		for(int l=0; l<this->p_ucell->atoms[it].nwl+1; l++)
		{
			for(int izeta=0; izeta<this->p_ucell->atoms[it].l_nchi[l]; izeta++)
			{
				double* ovlp_flzjlq_q = new double[GlobalV::NQX];
                double* qgrid = new double[GlobalV::NQX];
                for (int iq = 0; iq < GlobalV::NQX; iq++)
                {
                    qgrid[iq] = iq*GlobalV::DQ;
                }
                this->sbt.direct(l, 
								 this->n_rgrid[it][ic],
								 this->rgrid[it][ic].data(),
								 this->flz[it][ic].data(),
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
		for(int it = 0; it < this->p_ucell->ntype; it++)
		{
			std::stringstream ss;
			ss<<GlobalV::global_out_dir<<this->p_ucell->atoms[it].label<< "/LOCAL_G.dat";
			std::ofstream ofs(ss.str().c_str());
			for(int iq = 0; iq < GlobalV::NQX; iq++)
			{
				int ic=0;
				double energy_q = pow((double)iq*GlobalV::DQ, 2);
				ofs<<energy_q*this->p_ucell->tpiba2;
				for(int l = 0; l<this->p_ucell->atoms[it].nwl + 1; l++)
				{
					for(int N=0; N<this->p_ucell->atoms[it].l_nchi[l]; N++)
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

template <typename T, typename Device>
psi::Psi<T, Device>* psi_initializer_nao<T, Device>::cal_psig(int ik)
{
	ModuleBase::timer::tick("psi_initializer_nao", "initialize");
	assert(ik>=0);
	this->psig->fix_k(ik);
	const int npw = this->pw_wfc->npwk[ik];
	const int total_lm = ( this->p_ucell->lmax + 1) * ( this->p_ucell->lmax + 1);
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
	for (int it = 0; it < this->p_ucell->ntype; it++)
	{
		for (int ia = 0; ia < this->p_ucell->atoms[it].na; ia++)
		{
/* HERE LOOP OVER ALL ATOMIS */
            std::complex<double>* sk = this->sf->get_sk(ik, it, ia, this->pw_wfc);
            int ic = 0; // ic is a flatten index of chi, therefore it is defined here.
            for(int L = 0; L < this->p_ucell->atoms[it].nwl+1; L++)
			{
				std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, L); //mohan 2010-04-19
				for(int N=0; N < this->p_ucell->atoms[it].l_nchi[L]; N++)
				{
/* HERE LOOP OVER ALL NAOS */
					for(int ig=0; ig<npw; ig++)
					{
						ovlp_flzjlg[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(this->ovlp_flzjlq,
						it, ic, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * this->p_ucell->tpiba );
					}
/* FOR EVERY NAO IN EACH ATOM */
					if(GlobalV::NSPIN==4)
					{
/* FOR EACH SPIN CHANNEL */
						for(int is_N = 0; is_N < 2; is_N++) // rotate base
						//for(int is_N = 0; is_N < 1; is_N++)
						{
							if(L==0 && is_N==1) continue;
							if(this->p_ucell->atoms[it].ncpp.has_so)
							{
								const double j = std::abs(double(L+is_N) - 0.5);
								if (!(GlobalV::DOMAG||GlobalV::DOMAG_Z))
								{//atomic_wfc_so
									for(int m=0; m<2*L+1; m++)
									{
										const int lm = L*L+m;
										for(int ig=0; ig<npw; ig++)
										{
											//if(is_N==0)
											(*(this->psig))(ibasis, ig) = 
												this->template cast_to_T<T>(
													lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig]
												);
											(*(this->psig))(ibasis, ig + this->pw_wfc->npwk_max) =
												this->template cast_to_T<T>(
													lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig]
												);
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
												it, ic, GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * this->p_ucell->tpiba );

											chiaux[ig] += ovlp_flzjlg[ig] * (L+1.0) ;
											chiaux[ig] *= 1/(2.0*L+1.0);
										}
									}
									alpha = this->p_ucell->atoms[it].angle1[ia];
									gamma = -1 * this->p_ucell->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
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
											(*(this->psig))(ibasis, ig) = 
												this->template cast_to_T<T>(
													(cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup
												);
											(*(this->psig))(ibasis, ig + this->pw_wfc->npwk_max) =
												this->template cast_to_T<T>(
													(cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown
												);
                                            // second rotation with angle gamma around(OZ)
                                            fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                            fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI))*aux[ig];
											(*(this->psig))(ibasis+2*L+1,ig) =
												this->template cast_to_T<T>(
													(cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup
												);
											(*(this->psig))(ibasis+2*L+1, ig + this->pw_wfc->npwk_max) =
												this->template cast_to_T<T>(
													(cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown
												);
										}
                                        ibasis++;
                                    }
									ibasis += 2*L +1;
								} // end else INPUT.starting_spin_angle || !GlobalV::DOMAG
							} // end if this->p_ucell->atoms[it].has_so
							else
							{//atomic_wfc_nc
								double alpha, gamman;
								std::complex<double> fup, fdown;
								//alpha = this->p_ucell->magnet.angle1_[it];
								//gamman = -this->p_ucell->magnet.angle2_[it] + 0.5*ModuleBase::PI;
								alpha = this->p_ucell->atoms[it].angle1[ia];
								gamman = -this->p_ucell->atoms[it].angle2[ia] + 0.5*ModuleBase::PI;
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
										(*(this->psig))(ibasis, ig) = 
											this->template cast_to_T<T>(
												(cos(0.5*gamman) + ModuleBase::IMAG_UNIT*sin(0.5*gamman)) * fup
											);
										(*(this->psig))(ibasis, ig + this->pw_wfc->npwk_max) =
											this->template cast_to_T<T>(
												(cos(0.5*gamman) - ModuleBase::IMAG_UNIT*sin(0.5*gamman)) * fdown
											);
                                        // second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
										(*(this->psig))(ibasis+2*L+1,ig) =
											this->template cast_to_T<T>(
												(cos(0.5*gamman) + ModuleBase::IMAG_UNIT*sin(0.5*gamman)) * fup
											);
										(*(this->psig))(ibasis+2*L+1, ig + this->pw_wfc->npwk_max) =
											this->template cast_to_T<T>(
												(cos(0.5*gamman) - ModuleBase::IMAG_UNIT*sin(0.5*gamman)) * fdown
											);
                                    } // end ig
                                    ibasis++;
                                } // end m
								ibasis += 2*L+1;
							} // end else this->p_ucell->atoms[it].has_so
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
									this->template cast_to_T<T>(
										lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig]
									);
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
		this->random_t(this->psig->get_pointer(), ibasis, this->psig->get_nbands(), ik);
	}
	ModuleBase::timer::tick("psi_initializer_nao", "initialize");
	
	return this->psig;
}

template class psi_initializer_nao<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer_nao<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer_nao<double, psi::DEVICE_CPU>;
template class psi_initializer_nao<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer_nao<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer_nao<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer_nao<double, psi::DEVICE_GPU>;
template class psi_initializer_nao<float, psi::DEVICE_GPU>;
#endif
