#include "psi_initializer_atomic.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
// numerical algorithm support
#include "module_base/math_integral.h" // for numerical integration
#include "module_base/math_polyint.h" // for polynomial interpolation
#include "module_base/math_ylmreal.h" // for real spherical harmonics
#include "module_base/math_sphbes.h" // for spherical bessel functions
// basic functions support
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
// global variables definition
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
// io support
#include "module_io/write_pao.h"

// free function, compared with common radial function normalization, it does not multiply r to function
// due to pswfc is already multiplied by r
// template <typename T>
// void normalize(int n_rgrid, std::vector<T>& pswfcr, double* rab)
// {
//     std::vector<T> pswfc2r2(pswfcr.size());
//     std::transform(pswfcr.begin(), pswfcr.end(), pswfc2r2.begin(), [](T pswfc) { return pswfc * pswfc; });
//     T norm = ModuleBase::Integral::simpson(n_rgrid, pswfc2r2.data(), rab);
//     norm = sqrt(norm);
//     std::transform(pswfcr.begin(), pswfcr.end(), pswfcr.begin(), [norm](T pswfc) { return pswfc / norm; });
// }

template <typename T>
void psi_initializer_atomic<T>::allocate_ps_table()
{
   // find correct dimension for ovlp_flzjlq
    int dim1 = this->p_ucell_->ntype;
    int dim2 = 0; // dim2 should be the maximum number of pseudo atomic orbitals
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        dim2 = std::max(dim2, this->p_ucell_->atoms[it].ncpp.nchi);
    }
    if (dim2 == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_atomic<T>::allocate_table", "there is not ANY pseudo atomic orbital read in present system, recommand other methods, quit.");
    }
    int dim3 = PARAM.globalv.nqx;
    // allocate memory for ovlp_flzjlq
    this->ovlp_pswfcjlq_.create(dim1, dim2, dim3);
    this->ovlp_pswfcjlq_.zero_out();
}

template <typename T>
void psi_initializer_atomic<T>::initialize(const Structure_Factor* sf,         //< structure factor
                                           const ModulePW::PW_Basis_K* pw_wfc, //< planewave basis
                                           const UnitCell* p_ucell,            //< unit cell
                                           const K_Vectors* p_kv_in,
                                           const int& random_seed,       //< random seed
                                           const pseudopot_cell_vnl* p_pspot_nl,
                                           const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize");
    if(p_pspot_nl == nullptr)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_atomic<T>::initialize", 
                                 "pseudopot_cell_vnl object cannot be nullptr for atomic, quit.");
    }
    // import
    psi_initializer<T>::initialize(sf, pw_wfc, p_ucell, p_kv_in, random_seed, p_pspot_nl, rank);
    this->nbands_start_ = std::max(this->p_ucell_->natomwfc, PARAM.inp.nbands);
    this->nbands_complem_ = this->nbands_start_ - this->p_ucell_->natomwfc;
    // allocate
    this->allocate_ps_table();
    // then for generate random number to fill in the wavefunction
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize_only_once");
}

template <typename T>
void psi_initializer_atomic<T>::tabulate()
{
    if (PARAM.inp.use_paw)
    {
        return;
    }
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_ovlp_pswfcjlq");
    
    GlobalV::ofs_running << "\n Make real space PAO into reciprocal space." << std::endl;
    ModuleIO::print_PAOs(*this->p_ucell_);

    // Find the type of atom that has most mesh points.
    int max_msh = 0;
    for (int it=0; it<this->p_ucell_->ntype; it++)
    {
        max_msh = (this->p_ucell_->atoms[it].ncpp.msh > max_msh) ? this->p_ucell_->atoms[it].ncpp.msh : max_msh;
    }
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max mesh points in Pseudopotential",max_msh);
    
    this->ovlp_pswfcjlq_.zero_out();
    const int startq = 0;
    const double pref = ModuleBase::FOUR_PI / sqrt(this->p_ucell_->omega);
    std::vector<double> aux(max_msh);
    std::vector<double> vchi(max_msh);

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dq(describe PAO in reciprocal space)",PARAM.globalv.dq);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max q",PARAM.globalv.nqx);

    for (int it=0; it<this->p_ucell_->ntype; it++)
    {
		Atom* atom = &this->p_ucell_->atoms[it];

		GlobalV::ofs_running<<"\n number of pseudo atomic orbitals for "<<atom->label<<" is "<< atom->ncpp.nchi << std::endl;

        // QE uses atom->ncpp.mesh
        const int n_rgrid = (PARAM.inp.pseudo_mesh) ? atom->ncpp.mesh : atom->ncpp.msh;
        std::vector<double> chi2(n_rgrid);

        for (int ic = 0; ic < atom->ncpp.nchi ;ic++)
        {
            // check the unit condition
            for(int ir=0; ir<n_rgrid; ir++)
            {
                double chi = atom->ncpp.chi(ic, ir);
                chi2[ir] = chi * chi;
            }
            double unit = 0.0;
            ModuleBase::Integral::Simpson_Integral(n_rgrid, chi2.data(), atom->ncpp.rab.data(), unit);
            // liuyu add 2023-10-06
            if (unit < 1e-8)
            {
                // set occupancy to a small negative number so that this wfc
                // is not going to be used for starting wavefunctions
                atom->ncpp.oc[ic] = -1e-8;
                GlobalV::ofs_running << "WARNING: norm of atomic wavefunction # " << ic + 1 << " of atomic type "
                                     << atom->ncpp.psd << " is zero" << std::endl;
            }

            // only occupied states are normalized
            if (atom->ncpp.oc[ic] < 0)
            {
                continue;
            }

            // the US part if needed
            if (atom->ncpp.tvanp)
            {
                int kkbeta = atom->ncpp.kkbeta;
                if ((kkbeta % 2 == 0) && kkbeta > 0)
                {
                    kkbeta--;
                }
                std::vector<double> norm_beta(kkbeta);
                std::vector<double> work(atom->ncpp.nbeta);
                for (int ib = 0; ib < atom->ncpp.nbeta; ib++)
                {
                    bool match = false;
                    if (atom->ncpp.lchi[ic] == atom->ncpp.lll[ib])
                    {
                        if (atom->ncpp.has_so)
                        {
                            if (std::abs(atom->ncpp.jchi[ic] - atom->ncpp.jjj[ib]) < 1e-6)
                            {
                                match = true;
                            }
                        }
                        else
                        {
                            match = true;
                        }
                    }
                    if (match)
                    {
                        for (int ik = 0; ik < kkbeta; ik++)
                        {
                            norm_beta[ik] = atom->ncpp.betar(ib, ik) * atom->ncpp.chi(ic, ik);
                        }
                        ModuleBase::Integral::Simpson_Integral(kkbeta, norm_beta.data(), atom->ncpp.rab.data(), work[ib]);
                    }
                    else
                    {
                        work[ib] = 0.0;
                    }
                }
                for (int ib1 = 0; ib1 < atom->ncpp.nbeta; ib1++)
                {
                    for (int ib2 = 0; ib2 < atom->ncpp.nbeta; ib2++)
                    {
                        unit += atom->ncpp.qqq(ib1, ib2) * work[ib1] * work[ib2];
                    }
                }
            } // endif tvanp

            //=================================
            // normalize radial wave functions
            //=================================
            unit = std::sqrt(unit);
            if (std::abs(unit - 1.0) > 1e-6)
            {
                GlobalV::ofs_running << "WARNING: norm of atomic wavefunction # " << ic + 1 << " of atomic type "
                                     << atom->ncpp.psd << " is " << unit << ", renormalized" << std::endl;
                for (int ir = 0; ir < n_rgrid; ir++)
                {
                    atom->ncpp.chi(ic, ir) /= unit;
                }
            }

            const int l = atom->ncpp.lchi[ic];
            for (int iq = startq; iq < PARAM.globalv.nqx; iq++)
            {
                const double q = PARAM.globalv.dq * iq;
                ModuleBase::Sphbes::Spherical_Bessel(atom->ncpp.msh, atom->ncpp.r.data(), q, l, aux.data());
                for (int ir = 0; ir < atom->ncpp.msh; ir++)
                {
                    vchi[ir] = atom->ncpp.chi(ic, ir) * aux[ir] * atom->ncpp.r[ir];
                }

                double vqint = 0.0;
                ModuleBase::Integral::Simpson_Integral(atom->ncpp.msh, vchi.data(), atom->ncpp.rab.data(), vqint);

                this->ovlp_pswfcjlq_(it, ic, iq) = vqint * pref;
            }
        }
    }
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_ovlp_pswfcjlq");
}

std::complex<double> phase_factor(double arg, int mode)
{
    if(mode == 1) { return std::complex<double>(cos(arg),0); }
    else if (mode == -1) { return std::complex<double>(0, sin(arg)); }
    else if (mode == 0) { return std::complex<double>(cos(arg), sin(arg)); }
    else { return std::complex<double>(1,0); }
}

template <typename T>
void psi_initializer_atomic<T>::init_psig(T* psig,  const int& ik)
{
    ModuleBase::timer::tick("psi_initializer_atomic", "init_psig");
    const int npw = this->pw_wfc_->npwk[ik];
    const int npwk_max = this->pw_wfc_->npwk_max;
    int lmax = this->p_ucell_->lmax_ppwf;
    const int total_lm = (lmax + 1) * (lmax + 1);
    ModuleBase::matrix ylm(total_lm, npw);
    ModuleBase::GlobalFunc::ZEROS(psig, PARAM.globalv.npol * this->nbands_start_ * npwk_max);

    std::vector<std::complex<double>> aux(npw);
    std::vector<double> chiaux(npw);
    std::vector<ModuleBase::Vector3<double>> gk(npw);
    // I plan to use std::transform to replace the following for loop
    // but seems it is not as easy as I thought, the lambda function is not easy to write
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = this->pw_wfc_->getgpluskcar(ik, ig);
    }
    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk.data(), ylm);
    int index = 0;
    std::vector<double> ovlp_pswfcjlg(npw);
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        for (int ia = 0; ia < this->p_ucell_->atoms[it].na; ia++)
        {
/* FOR EVERY ATOM */
            // I think it is always a BAD idea to new one pointer in a function, then return it
            // it indicates the ownership of the pointer and behind memory is transferred to the caller
            // then one must manually delete it, makes new-delete not symmetric
            std::complex<double> *sk = this->sf_->get_sk(ik, it, ia, this->pw_wfc_);
            for (int ipswfc = 0; ipswfc < this->p_ucell_->atoms[it].ncpp.nchi; ipswfc++)
            {
/* FOR EVERY PSWFC OF ATOM */
                if (this->p_ucell_->atoms[it].ncpp.oc[ipswfc] >= 0.0)
                {
/* IF IS OCCUPIED, GET L */
                    const int l = this->p_ucell_->atoms[it].ncpp.lchi[ipswfc];
                    std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, l);

                    for (int ig=0; ig<npw; ig++)
                    {
                        ovlp_pswfcjlg[ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
                            this->ovlp_pswfcjlq_, it, ipswfc, 
                            PARAM.globalv.nqx, PARAM.globalv.dq, gk[ig].norm() * this->p_ucell_->tpiba );
                    }
/* NSPIN == 4 */
                    if(PARAM.inp.nspin == 4)
                    {
                        if(this->p_ucell_->atoms[it].ncpp.has_so)
                        {
                            Soc soc; soc.rot_ylm(l + 1);
                            const double j = this->p_ucell_->atoms[it].ncpp.jchi[ipswfc];
    /* NOT NONCOLINEAR CASE, rotation matrix become identity */
                            if (!(PARAM.globalv.domag||PARAM.globalv.domag_z))
                            {
                                double cg_coeffs[2];
                                for(int m = -l-1; m < l+1; m++)
                                {
                                    cg_coeffs[0] = soc.spinor(l, j, m, 0);
                                    cg_coeffs[1] = soc.spinor(l, j, m, 1);
                                    if (fabs(cg_coeffs[0]) > 1e-8 || fabs(cg_coeffs[1]) > 1e-8)
                                    {
                                        for(int is = 0; is < 2; is++)
                                        {
                                            if(fabs(cg_coeffs[is]) > 1e-8)
                                            {
        /* GET COMPLEX SPHERICAL HARMONIC FUNCTION */
                                                const int ind = this->p_pspot_nl_->lmaxkb + soc.sph_ind(l,j,m,is); // ind can be l+m, l+m+1, l+m-1
                                                std::fill(aux.begin(), aux.end(), std::complex<double>(0.0, 0.0));
                                                for(int n1 = 0; n1 < 2*l+1; n1++)
                                                {
                                                    const int lm = l*l +n1;
                                                    std::complex<double> umM = soc.rotylm(n1, ind);
                                                    if(std::abs(umM) > 1e-8)
                                                    {
                                                        for(int ig = 0; ig < npw; ig++)
                                                        {
                                                            aux[ig] += umM * ylm(lm, ig);
                                                        }
                                                    }
                                                }
                                                for(int ig = 0; ig < npw; ig++)
                                                {
                                                    psig[(2 * index + is) * npwk_max + ig] = this->template cast_to_T<T>(
                                                        lphase * cg_coeffs[is] * sk[ig] * aux[ig] * ovlp_pswfcjlg[ig]);
                                                }
                                            }
                                            else
                                            {
                                                for (int ig = 0; ig < npw; ig++)
                                                {
                                                    psig[(2 * index + is) * npwk_max + ig]
                                                        = this->template cast_to_T<T>(std::complex<double>(0.0, 0.0));
                                                }
                                            }
                                        }
                                        index++;
                                    }
                                }
                            }
                            else
                            {
    /* NONCONLINEAR CASE, will use [[cos(a/2)*exp(-ib/2), sin(a/2)*exp(ib/2)], [-sin(a/2)*exp(-ib/2), cos(a/2)*exp(ib/2)]] to rotate */
                                int ipswfc_noncolin_soc=0;
        /* J = L - 1/2 -> continue */
        /* J = L + 1/2 */
								if(fabs(j - l + 0.5) < 1e-4) 
								{
									continue;
								}
								chiaux.clear(); 
								chiaux.resize(npw);
        /* L == 0 */
								if(l == 0) 
								{
									std::memcpy(chiaux.data(), ovlp_pswfcjlg.data(), npw * sizeof(double));
								}
                                else
                                {
        /* L != 0, scan pswfcs that have the same L and satisfy J(pswfc) = L - 0.5 */
                                    for(int jpsiwfc = 0; jpsiwfc < this->p_ucell_->atoms[it].ncpp.nchi; jpsiwfc++)
                                    {
                                        if(
                                            (this->p_ucell_->atoms[it].ncpp.lchi[jpsiwfc] == l)
                                          &&(fabs(this->p_ucell_->atoms[it].ncpp.jchi[jpsiwfc] - l + 0.5) < 1e-4))
                                        {
                                            ipswfc_noncolin_soc = jpsiwfc;
                                            break;
                                        }
                                    }
                                    for(int ig=0;ig<npw;ig++)
                                    {
            /* average <pswfc_a|jl(q)> and <pswfc_b(j=l-1/2)|jl(q)>, a and b seem not necessarily to be equal */
                                        chiaux[ig] =  l *
                                            ModuleBase::PolyInt::Polynomial_Interpolation(
                                                this->ovlp_pswfcjlq_, it, ipswfc_noncolin_soc, 
                                                PARAM.globalv.nqx, PARAM.globalv.dq, gk[ig].norm() * this->p_ucell_->tpiba);
                                        chiaux[ig] += ovlp_pswfcjlg[ig] * (l + 1.0) ;
                                        chiaux[ig] *= 1/(2.0*l+1.0);
                                    }
                                }
            /* ROTATE ACCORDING TO NONCOLINEAR */
                                double alpha = this->p_ucell_->atoms[it].angle1[ia];
                                double gamma = -1 * this->p_ucell_->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
                                std::complex<double> fup, fdw;

                                for(int m = 0; m < 2*l+1; m++)
                                {
                                    const int lm = l*l +m;
                                    if(index+2*l+1 > this->p_ucell_->natomwfc)
                                    {
                                        std::cout<<__FILE__<<__LINE__<<" "<<index<<" "<<this->p_ucell_->natomwfc<<std::endl;
                                        //ModuleBase::WARNING_QUIT("psi_initializer_atomic<T>::init_psig()","error: too many wfcs");
                                    }
                                    for(int ig = 0;ig<npw;ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm,ig) * chiaux[ig];
                                    }
                                    //rotate wfc as needed
                                    //first rotation with angle alpha around (OX)
                                    for(int ig = 0;ig<npw;ig++)
                                    {
                                        fup = phase_factor(0.5*alpha,  1)*aux[ig];
                                        fdw = phase_factor(0.5*alpha, -1)*aux[ig];
                                        //build the orthogonal wfc
                                        //first rotation with angle (alpha + ModuleBase::PI) around (OX)
                                        psig[index * 2 * npwk_max + ig]
                                            = this->template cast_to_T<T>(phase_factor(0.5 * gamma, 0) * fup);
                                        psig[(index * 2 + 1) * npwk_max + ig]
                                            = this->template cast_to_T<T>(phase_factor(-0.5 * gamma, 0) * fdw);
                                        //second rotation with angle gamma around(OZ)
                                        fup = phase_factor(0.5*(alpha + ModuleBase::PI),  1)*aux[ig];
                                        fdw = phase_factor(0.5*(alpha + ModuleBase::PI), -1)*aux[ig];
                                        psig[(index + 2 * l + 1) * 2 * npwk_max + ig]
                                            = this->template cast_to_T<T>(phase_factor(0.5 * gamma, 0) * fup);
                                        psig[((index + 2 * l + 1) * 2 + 1) * npwk_max + ig]
                                            = this->template cast_to_T<T>(phase_factor(-0.5 * gamma, 0) * fdw);
                                    }
                                    index++;
                                }
                                index += 2*l +1;
                            }
                        }
                        else
                        {//atomic_wfc_nc
                            double alpha=0.0;
                            double gamman=0.0;
                            std::complex<double> fup, fdown;
                            //alpha = this->p_ucell_->magnet.angle1_[it];
                            //gamman = -this->p_ucell_->magnet.angle2_[it] + 0.5*ModuleBase::PI;
                            alpha = this->p_ucell_->atoms[it].angle1[ia];
                            gamman = -1 * this->p_ucell_->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
                            for(int m = 0; m < 2*l+1; m++)
                            {
                                const int lm = l*l +m;
                                if(index+2*l+1 > this->p_ucell_->natomwfc)
                                {
                                    std::cout<<__FILE__<<__LINE__<<" "<<index<<" "<<this->p_ucell_->natomwfc<<std::endl;
                                    //ModuleBase::WARNING_QUIT("psi_initializer_atomic<T>::init_psig()","error: too many wfcs");
                                }
                                for(int ig = 0;ig<npw;ig++)
                                {
                                     aux[ig] = sk[ig] * ylm(lm,ig) * ovlp_pswfcjlg[ig];
                                }
                                //rotate function
                                //first, rotation with angle alpha around(OX)
                                for(int ig = 0; ig<npw; ig++)
                                {
                                    fup = cos(0.5 * alpha) * aux[ig];
                                    fdown = ModuleBase::IMAG_UNIT * sin(0.5 * alpha) * aux[ig];
                                    // build the orthogonal wfc
                                    // first rotation with angle(alpha+ModuleBase::PI) around(OX)
                                    psig[index * 2 * npwk_max + ig] = this->template cast_to_T<T>(
                                        (cos(0.5 * gamman) + ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fup);
                                    psig[(index * 2 + 1) * npwk_max + ig] = this->template cast_to_T<T>(
                                        (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown);
                                    // second rotation with angle gamma around(OZ)
                                    fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                    fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                    psig[(index + 2 * l + 1) * 2 * npwk_max + ig] = this->template cast_to_T<T>(
                                        (cos(0.5 * gamman) + ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fup);
                                    psig[((index + 2 * l + 1) * 2 + 1) * npwk_max + ig] = this->template cast_to_T<T>(
                                        (cos(0.5 * gamman) - ModuleBase::IMAG_UNIT * sin(0.5 * gamman)) * fdown);
                                }
                                index++;
                            }
                            index += 2*l+1;
                        }
                    }
                    else
                    {
                        for (int m = 0; m < 2*l+1; m++)
                        {
                            const int lm = l * l + m;
                            for (int ig = 0; ig < npw; ig++)
                            {
                                psig[index * npwk_max + ig]
                                    = this->template cast_to_T<T>(lphase * sk[ig] * ylm(lm, ig) * ovlp_pswfcjlg[ig]);
                            }
                            index++;
                        }
                    }
                }
            }
			delete [] sk;
        }
    }
	/* complement the rest of bands if there are */
	if(this->nbands_complem() > 0)
	{
		this->random_t(psig, index, this->nbands_start_, ik);
	}
    ModuleBase::timer::tick("psi_initializer_atomic", "init_psig");
}

template class psi_initializer_atomic<std::complex<double>>;
template class psi_initializer_atomic<std::complex<float>>;
// gamma point calculation
template class psi_initializer_atomic<double>;
template class psi_initializer_atomic<float>;
