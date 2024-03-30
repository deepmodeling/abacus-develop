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
#include "module_base/parallel_reduce.h"
#endif

/*
I don't know why some variables are distributed while others not... for example the orbital_files...
We need not only read and import, but also distribute here
*/

// free function, not needed to be a member of psi_initializer_nao
void normalize(const std::vector<double>& r, std::vector<double>& flz)
{
    std::vector<double> flz2r2(r.size());
    std::transform(r.begin(), r.end(), flz.begin(), flz2r2.begin(), [](double r, double flz){return flz*flz*r*r;});
    double dr = r[1] - r[0];
    double norm = ModuleBase::Integral::simpson(r.size(), flz2r2.data(), dr);
    norm = sqrt(norm);
    std::transform(flz.begin(), flz.end(), flz.begin(), [norm](double flz){return flz/norm;});
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::read_external_orbs(std::string* orbital_files,
                                                        const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_nao", "read_external_orbs");
    if(rank == 0)
    {
        for (int itype = 0; itype < this->p_ucell_->ntype; itype++)
        {
            this->orbital_files_.push_back(orbital_files[itype]);
        }
        for(int it = 0; it < this->p_ucell_->ntype; it++)
        {
            // number of chi per atomtype
            int nchi = 0;
            for(int l = 0; l <= this->p_ucell_->atoms[it].nwl; l++)
            {
                nchi += this->p_ucell_->atoms[it].l_nchi[l];
            }
            
            std::vector<int> n_rgrid_it;
            std::vector<std::vector<double>> rgrid_it;
            std::vector<std::vector<double>> rvalue_it;

            std::ifstream ifs_it;
            ifs_it.open(GlobalV::global_orbital_dir+this->orbital_files_[it]);

            if(!ifs_it)
            {
                GlobalV::ofs_warning<<"psi_initializer_nao<T, Device>::read_orbital_files: cannot open orbital file: "<<this->orbital_files_[it]<<std::endl;
                ModuleBase::WARNING_QUIT("psi_initializer_nao<T, Device>::read_orbital_files", "cannot open orbital file.");
            }
            else
            {
                GlobalV::ofs_running<<"psi_initializer_nao<T, Device>::read_orbital_files: reading orbital file: "<<this->orbital_files_[it]<<std::endl;
            }
            ifs_it.close();

            int ichi_overall = 0;
            // check nwl and nchi for each rank
            
            for(int l = 0; l <= this->p_ucell_->atoms[it].nwl; l++)
            {
                for(int ichi = 0; ichi < this->p_ucell_->atoms[it].l_nchi[l]; ichi++)
                {
                    int n_rgrid_ichi;
                    std::vector<double> rgrid_ichi;
                    std::vector<double> rvalue_ichi;

                    GlobalV::ofs_running<<"-------------------------------------- "<<std::endl;
                    GlobalV::ofs_running<<" reading orbital of element "<<this->p_ucell_->atoms[it].label<<std::endl
                                        <<" angular momentum l = "<<l<<std::endl
                                        <<" index of chi = "<<ichi<<std::endl;

                    ifs_it.open(GlobalV::global_orbital_dir+this->orbital_files_[it]);
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
                            GlobalV::ofs_warning<<" psi_initializer_nao<T, Device>::read_orbital_files: cannot find orbital of element "<<this->p_ucell_->atoms[it].label<<std::endl
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
                                double rvalue_ichi_ir;
                                ifs_it>>rvalue_ichi_ir;
                                rvalue_ichi.push_back(rvalue_ichi_ir);
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
                    // before push back, normalize the rvalue_ichi, 2024/03/19, kirk0830
                    // turn off normalize, 2024/03/22, kirk0830
                    //normalize(rgrid_ichi, rvalue_ichi);
                    rvalue_it.push_back(rvalue_ichi);
                    ++ichi_overall;
                }
            }
            this->n_rgrid_.push_back(n_rgrid_it);
            this->rgrid_.push_back(rgrid_it);
            this->rvalue_.push_back(rvalue_it);
            GlobalV::ofs_running<<"-------------------------------------- "<<std::endl;
        }
    }
    // MPI additional implementation
    #ifdef __MPI
    // bcast fname
    if(rank != 0) this->orbital_files_.resize(this->p_ucell_->ntype);
    Parallel_Common::bcast_string(this->orbital_files_.data(), this->p_ucell_->ntype);
    // bcast orbital data
    // resize
    if(rank != 0) this->n_rgrid_.resize(this->p_ucell_->ntype);
    int nchi[this->p_ucell_->ntype];
    if(rank == 0)
    {
        for(int it = 0; it < this->p_ucell_->ntype; it++) nchi[it] = this->n_rgrid_[it].size();
    }
    // bcast
    Parallel_Common::bcast_int(nchi, this->p_ucell_->ntype);
    // resize
    if(rank != 0)
    {
        this->n_rgrid_.resize(this->p_ucell_->ntype);
        this->rgrid_.resize(this->p_ucell_->ntype);
        this->rvalue_.resize(this->p_ucell_->ntype);
        for(int it = 0; it < this->p_ucell_->ntype; it++)
        {
            this->n_rgrid_[it].resize(nchi[it]);
            this->rgrid_[it].resize(nchi[it]);
            this->rvalue_[it].resize(nchi[it]);
        }
    }
    // bcast
    for(int it = 0; it < this->p_ucell_->ntype; it++) Parallel_Common::bcast_int(this->n_rgrid_[it].data(), nchi[it]);
    // resize
    if(rank != 0)
    {
        for(int it = 0; it < this->p_ucell_->ntype; it++)
        {
            for(int ichi = 0; ichi < nchi[it]; ichi++)
            {
                this->rgrid_[it][ichi].resize(this->n_rgrid_[it][ichi]);
                this->rvalue_[it][ichi].resize(this->n_rgrid_[it][ichi]);
            }
        }
    }
    // bcast
    for(int it = 0; it < this->p_ucell_->ntype; it++)
    {
        for(int ichi = 0; ichi < nchi[it]; ichi++)
        {
            Parallel_Common::bcast_double(this->rgrid_[it][ichi].data(), this->n_rgrid_[it][ichi]);
            Parallel_Common::bcast_double(this->rvalue_[it][ichi].data(), this->n_rgrid_[it][ichi]);
        }
    }
    #endif
    ModuleBase::timer::tick("psi_initializer_nao", "read_external_orbs");
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::allocate_table()
{
    // find correct dimension for ovlp_flzjlq
    int dim1 = this->p_ucell_->ntype;
    int dim2 = 0; // dim2 should be the maximum number of zeta for each atomtype
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        int nzeta = 0;
        for (int l = 0; l < this->p_ucell_->atoms[it].nwl+1; l++)
        {
            nzeta += this->p_ucell_->atoms[it].l_nchi[l];
        }
        dim2 = (nzeta > dim2) ? nzeta : dim2;
    }
    if (dim2 == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_nao<T, Device>::psi_initializer_nao", "there is not ANY numerical atomic orbital read in present system, quit.");
    }
    int dim3 = GlobalV::NQX;
    // allocate memory for ovlp_flzjlq
    this->ovlp_flzjlq_.create(dim1, dim2, dim3);
    this->ovlp_flzjlq_.zero_out();
}


#ifdef __MPI
template <typename T, typename Device>
void psi_initializer_nao<T, Device>::initialize(Structure_Factor* sf,
                                                ModulePW::PW_Basis_K* pw_wfc,
                                                UnitCell* p_ucell,
                                                Parallel_Kpoints* p_parakpts,
                                                const int& random_seed,
                                                pseudopot_cell_vnl* p_pspot_nl,
                                                const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_nao", "initialize_mpi");
    // import
    this->sf_ = sf;
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_parakpts_ = p_parakpts;
    this->p_pspot_nl_ = p_pspot_nl;
    this->random_seed_ = random_seed;
    // allocate
    this->allocate_table();
    this->read_external_orbs(this->p_ucell_->orbital_fn, rank);
    //this->cal_ovlp_flzjlq(); //because GlobalV::NQX will change during vcrelax, so it should be called in both init and init_after_vc
    ModuleBase::timer::tick("psi_initializer_nao", "initialize_mpi");
}
#else
template <typename T, typename Device>
void psi_initializer_nao<T, Device>::initialize(Structure_Factor* sf,
                                                ModulePW::PW_Basis_K* pw_wfc,
                                                UnitCell* p_ucell,
                                                const int& random_seed,
                                                pseudopot_cell_vnl* p_pspot_nl)
{
    ModuleBase::timer::tick("psi_initializer_nao", "initialize_serial");
    // import
    this->sf_ = sf;
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_pspot_nl_ = p_pspot_nl;
    this->random_seed_ = random_seed;
    // allocate
    this->allocate_table();
    this->read_external_orbs(this->p_ucell_->orbital_fn, 0);
    //this->cal_ovlp_flzjlq(); //because GlobalV::NQX will change during vcrelax, so it should be called in both init and init_after_vc
    ModuleBase::timer::tick("psi_initializer_nao", "initialize_serial");
}
#endif

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::tabulate()
{
    ModuleBase::timer::tick("psi_initializer_nao", "tabulate");
    this->ovlp_flzjlq_.zero_out();
    for(int it=0; it<this->p_ucell_->ntype; it++)
    {
        int ic=0;
        for(int l=0; l<this->p_ucell_->atoms[it].nwl+1; l++)
        {
            for(int izeta=0; izeta<this->p_ucell_->atoms[it].l_nchi[l]; izeta++)
            {
                std::vector<double> ovlp_flzjlq_q(GlobalV::NQX);
                std::vector<double> qgrid(GlobalV::NQX);
                for (int iq = 0; iq < GlobalV::NQX; iq++)
                {
                    qgrid[iq] = iq*GlobalV::DQ;
                }
                this->sbt.direct(l, 
                                 this->n_rgrid_[it][ic],
                                 this->rgrid_[it][ic].data(),
                                 this->rvalue_[it][ic].data(),
                                 GlobalV::NQX, 
                                 qgrid.data(), 
                                 ovlp_flzjlq_q.data());
                for(int iq = 0; iq < GlobalV::NQX; iq++)
                {
                    this->ovlp_flzjlq_(it, ic, iq) = ovlp_flzjlq_q[iq];
                }
                ++ic;
            }
        }
    }
    ModuleBase::timer::tick("psi_initializer_nao", "tabulate");
}

template <typename T, typename Device>
void psi_initializer_nao<T, Device>::proj_ao_onkG(int ik)
{
    ModuleBase::timer::tick("psi_initializer_nao", "initialize");
    assert(ik>=0);
    this->psig_->fix_k(ik);
    const int npw = this->pw_wfc_->npwk[ik];
    const int total_lm = ( this->p_ucell_->lmax + 1) * ( this->p_ucell_->lmax + 1);
    ModuleBase::matrix ylm(total_lm, npw);

    std::vector<std::complex<double>> aux(npw);
    std::vector<ModuleBase::Vector3<double>> gk(npw);
    for(int ig=0;ig<npw;ig++)
    {
        gk[ig] = this->pw_wfc_->getgpluskcar(ik, ig);
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk.data(), ylm);
    //int index = 0;
    std::vector<double> ovlp_flzjlg(npw);
    int ibasis=0;
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
/* HERE LOOP OVER ALL TYPES */
        for (int ia = 0; ia < this->p_ucell_->atoms[it].na; ia++)
        {
/* HERE LOOP OVER ALL ATOMS */
            std::complex<double>* sk = this->sf_->get_sk(ik, it, ia, this->pw_wfc_);
            int ic = 0; // ic is a flatten index of chi, therefore it is defined here.
            for(int L = 0; L < this->p_ucell_->atoms[it].nwl+1; L++)
            {
                std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, L); //mohan 2010-04-19
                for(int N=0; N < this->p_ucell_->atoms[it].l_nchi[L]; N++)
                {
/* HERE LOOP OVER ALL NAOS */
					/* 
						for already using flattened 1d index of chi, which folds l and n, the spherical bessel
						transformation of numerical orbital function, is indiced by it and ic, is needed to
						interpolate everytime when ic updates, therefore everytime when present orbital is done
					*/
                    for(int ig=0; ig<npw; ig++)
                    {
                        ovlp_flzjlg[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
							this->ovlp_flzjlq_, // the spherical bessel transform of numerical orbital function
                        	it, ic, 		   // each (it, ic)-pair defines a unique numerical orbital function
							GlobalV::NQX, GlobalV::DQ, // grid number and grid spacing of q
							gk[ig].norm() * this->p_ucell_->tpiba // norm of (G+k) = K
							);
                    }
/* FOR EVERY NAO IN EACH ATOM */
                    if(GlobalV::NSPIN==4)
                    {
/* FOR EACH SPIN CHANNEL */
                        for(int is_N = 0; is_N < 2; is_N++) // rotate base
                        //for(int is_N = 0; is_N < 1; is_N++)
                        {
                            if(L==0 && is_N==1)
							{
								continue;
							}
							else
							{
                                const double j = fabs(double(L+is_N) - 0.5);
                                double alpha, gamma;
                                std::complex<double> fup,fdown;
                                if(fabs(j - L + 0.5) < 1e-4)
                                {
                                    continue;
                                }
                                alpha = this->p_ucell_->atoms[it].angle1[ia];
                                gamma = -1 * this->p_ucell_->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
                                for(int m = 0;m<2*L+1;m++)
                                {
                                    const int lm = L*L + m;
                                    for (int ig = 0; ig < npw; ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm,ig) * ovlp_flzjlg[ig];
                                    }
                                    for(int ig = 0;ig<npw;ig++)
                                    {
                                        fup = cos(0.5 * alpha) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                        //build the orthogonal wfc
                                        //first rotation with angle (alpha + ModuleBase::PI) around (OX)
                                        (*(this->psig_))(ibasis, ig) = 
                                            this->template cast_to_T<T>(
                                                (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup
                                            );
                                        (*(this->psig_))(ibasis, ig + this->pw_wfc_->npwk_max) =
                                            this->template cast_to_T<T>(
                                                (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown
                                            );
                                        // second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI))*aux[ig];
                                        (*(this->psig_))(ibasis+2*L+1,ig) =
                                            this->template cast_to_T<T>(
                                                (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup
                                            );
                                       (*(this->psig_))(ibasis+2*L+1, ig + this->pw_wfc_->npwk_max) =
                                            this->template cast_to_T<T>(
                                                (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown
                                            );
                                    }
                                    ibasis++;
                                }
                                ibasis += 2*L +1;
							}
                        } // end for is_N
                    } // end if GlobalV::NONCOLIN
                    else{//LSDA and nomagnet case
    /* DOES NOT DISTINGUISH m QUANTUM NUMBER FOR CHI */
                        for(int m = 0; m < 2*L+1; m++)
                        {
                            const int lm = L*L+m;
                            for(int ig=0; ig<npw; ig++)
                            {
                                (*(this->psig_))(ibasis, ig) =  this->template cast_to_T<T>(lphase * sk[ig] * ylm(lm, ig) * ovlp_flzjlg[ig]);
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
    /* complement the rest of bands if there are */
    if(this->nbands_complem() > 0)
    {
        this->random_t(this->psig_->get_pointer(), ibasis, this->psig_->get_nbands(), ik);
    }
    ModuleBase::timer::tick("psi_initializer_nao", "initialize");
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
