#include "psi_initializer_atomic.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
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

// free function, compared with common radial function normalization, it does not multiply r to function
// due to pswfc is already multiplied by r
template <typename T>
void normalize(int n_rgrid, std::vector<T>& pswfcr, double* rab)
{
    std::vector<T> pswfc2r2(pswfcr.size());
    std::transform(pswfcr.begin(), pswfcr.end(), pswfc2r2.begin(), [](T pswfc) { return pswfc * pswfc; });
    T norm = ModuleBase::Integral::simpson(n_rgrid, pswfc2r2.data(), rab);
    norm = sqrt(norm);
    std::transform(pswfcr.begin(), pswfcr.end(), pswfcr.begin(), [norm](T pswfc) { return pswfc / norm; });
}

template <typename T, typename Device>
void psi_initializer_atomic<T, Device>::allocate_table()
{
   // find correct dimension for ovlp_flzjlq
    int dim1 = this->p_ucell_->ntype;
    int dim2 = 0; // dim2 should be the maximum number of pseudo atomic orbitals
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        dim2 = (this->p_ucell_->atoms[it].ncpp.nchi > dim2) ? this->p_ucell_->atoms[it].ncpp.nchi : dim2;
    }
    if (dim2 == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_atomic<T, Device>::allocate_table", "there is not ANY pseudo atomic orbital read in present system, recommand other methods, quit.");
    }
    int dim3 = GlobalV::NQX;
    // allocate memory for ovlp_flzjlq
    this->ovlp_pswfcjlq_.create(dim1, dim2, dim3);
    this->ovlp_pswfcjlq_.zero_out();
}

#ifdef __MPI
template <typename T, typename Device>
void psi_initializer_atomic<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                   ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                   UnitCell* p_ucell,                                   //< unit cell
                                                   Parallel_Kpoints* p_parakpts,                        //< parallel kpoints
                                                   const int& random_seed,                          //< random seed
                                                   pseudopot_cell_vnl* p_pspot_nl,
                                                   const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize");
    if(p_pspot_nl == nullptr) ModuleBase::WARNING_QUIT("psi_initializer_atomic<T, Device>::initialize_only_once", 
                                                       "pseudopot_cell_vnl object cannot be mullptr for atomic, quit.");
    // import
    this->sf_ = sf;
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_parakpts_ = p_parakpts;
    this->p_pspot_nl_ = p_pspot_nl;
    this->random_seed_ = random_seed;
    // allocate
    this->allocate_table();
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize_only_once");
}
#else
template <typename T, typename Device>
void psi_initializer_atomic<T, Device>::initialize(Structure_Factor* sf,                                //< structure factor
                                                   ModulePW::PW_Basis_K* pw_wfc,                        //< planewave basis
                                                   UnitCell* p_ucell,                                   //< unit cell
                                                   const int& random_seed,                          //< random seed
                                                   pseudopot_cell_vnl* p_pspot_nl)
{
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize");
    if(p_pspot_nl == nullptr) ModuleBase::WARNING_QUIT("psi_initializer_atomic<T, Device>::initialize_only_once", 
                                                       "pseudopot_cell_vnl object cannot be mullptr for atomic, quit.");
    // import
    this->sf_ = sf;
    this->pw_wfc_ = pw_wfc;
    this->p_ucell_ = p_ucell;
    this->p_pspot_nl_ = p_pspot_nl;
    this->random_seed_ = random_seed;
    // allocate
    this->allocate_table();
    ModuleBase::timer::tick("psi_initializer_atomic", "initialize_only_once");
}
#endif

template <typename T, typename Device>
void psi_initializer_atomic<T, Device>::tabulate()
{
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_ovlp_pswfcjlq");
    int maxn_rgrid = 0;
    double* qgrid = new double[GlobalV::NQX];
    for (int iq = 0; iq < GlobalV::NQX; iq++)
    {
        qgrid[iq] = GlobalV::DQ * iq;
    }
    for (int it=0; it<this->p_ucell_->ntype; it++)
    {
        maxn_rgrid = (this->p_ucell_->atoms[it].ncpp.msh > maxn_rgrid) ? this->p_ucell_->atoms[it].ncpp.msh : maxn_rgrid;
    }
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max mesh points in Pseudopotential",maxn_rgrid);
    
    const double pref = ModuleBase::FOUR_PI / sqrt(this->p_ucell_->omega);

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dq(describe PAO in reciprocal space)",GlobalV::DQ);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max q",GlobalV::NQX);

    for (int it=0; it<this->p_ucell_->ntype; it++)
    {
		Atom* atom = &this->p_ucell_->atoms[it];

		GlobalV::ofs_running<<"\n number of pseudo atomic orbitals for "<<atom->label<<" is "<< atom->ncpp.nchi << std::endl;

        for (int ic = 0; ic < atom->ncpp.nchi ;ic++)
        {
            int n_rgrid = (GlobalV::PSEUDO_MESH)?atom->ncpp.mesh:atom->ncpp.msh;
            std::vector<double> pswfcr(n_rgrid);
            for (int ir=0; ir<n_rgrid; ir++) pswfcr[ir] = atom->ncpp.chi(ic, ir);
            normalize(n_rgrid, pswfcr, atom->ncpp.rab);
            if (atom->ncpp.oc[ic] >= 0.0) // reasonable occupation number, but is it always true?
            {
                const int l = atom->ncpp.lchi[ic];
                double* ovlp_pswfcjlq_q = new double[GlobalV::NQX];
                this->sbt.direct(l, atom->ncpp.msh, atom->ncpp.r, pswfcr.data(), GlobalV::NQX, qgrid, ovlp_pswfcjlq_q, 1);
                for (int iq = 0; iq < GlobalV::NQX; iq++)
                {
                    this->ovlp_pswfcjlq_(it, ic, iq) = pref * ovlp_pswfcjlq_q[iq];
                }
                delete [] ovlp_pswfcjlq_q;
            }
        }
    }
    delete [] qgrid;
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_ovlp_pswfcjlq");
}

std::complex<double> phase_factor(double arg, int mode)
{
    if(mode == 1) return std::complex<double>(cos(arg),0);
    else if (mode == -1) return std::complex<double>(0, sin(arg));
    else if (mode == 0) return std::complex<double>(cos(arg), sin(arg));
    else return std::complex<double>(1,0);
}

template <typename T, typename Device>
psi::Psi<T, Device>* psi_initializer_atomic<T, Device>::cal_psig(int ik)
{
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_psig");
    this->psig->fix_k(ik);
    //this->print_status(psi);
    const int npw = this->pw_wfc_->npwk[ik];
    int lmax = this->p_ucell_->lmax_ppwf;
    const int total_lm = (lmax + 1) * (lmax + 1);
    ModuleBase::matrix ylm(total_lm, npw);
    std::complex<double> *aux = new std::complex<double>[npw];
    double *chiaux = nullptr;

    ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3 <double> [npw];
    for (int ig = 0; ig < npw; ig++)
    {
        gk[ig] = this->pw_wfc_->getgpluskcar(ik, ig);
    }
    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk, ylm);
    int index = 0;
    double *ovlp_pswfcjlg = new double[npw];
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        for (int ia = 0; ia < this->p_ucell_->atoms[it].na; ia++)
        {
/* FOR EVERY ATOM */
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
                            GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * this->p_ucell_->tpiba );
                    }
/* NSPIN == 4 */
                    if(GlobalV::NSPIN == 4)
                    {
                        if(this->p_ucell_->atoms[it].ncpp.has_so)
                        {
                            Soc soc; soc.rot_ylm(l + 1);
                            const double j = this->p_ucell_->atoms[it].ncpp.jchi[ipswfc];
    /* NOT NONCOLINEAR CASE, rotation matrix become identity */
                            if (!(GlobalV::DOMAG||GlobalV::DOMAG_Z))
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
                                                ModuleBase::GlobalFunc::ZEROS(aux, npw);
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
                                                    (*(this->psig))(index, 
                                                                    ig + this->pw_wfc_->npwk_max*is ) =
                                                                        this->template cast_to_T<T>(
                                                                            lphase * cg_coeffs[is] * sk[ig] * aux[ig] * ovlp_pswfcjlg[ig]
                                                                        );
                                                }
                                            }
                                            else
                                            {
                                                for(int ig=0; ig < npw; ig++)
                                                {
                                                    (*(this->psig))(index, 
                                                                    ig + this->pw_wfc_->npwk_max*is ) = 
                                                                        this->template cast_to_T<T>(
                                                                            std::complex<double>(0.0, 0.0)
                                                                            );
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
                                int ipswfc_noncolin_soc;
        /* J = L - 1/2 -> continue */
        /* J = L + 1/2 */
                                if(fabs(j - l + 0.5) < 1e-4) continue;
                                delete[] chiaux; chiaux = new double [npw];
        /* L == 0 */
                                if(l == 0)
                                {
                                    for(int ig = 0; ig < npw; ig++)
                                    {
                                        chiaux[ig] = ovlp_pswfcjlg[ig];
                                    }
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
                                                GlobalV::NQX, GlobalV::DQ, gk[ig].norm() * this->p_ucell_->tpiba);
                                        chiaux[ig] += ovlp_pswfcjlg[ig] * (l + 1.0) ;
                                        chiaux[ig] *= 1/(2.0*l+1.0);
                                    }
                                }
            /* ROTATE ACCORDING TO NONCOLINEAR */
                                double alpha, gamma;
                                std::complex<double> fup, fdw;
                                alpha = this->p_ucell_->atoms[it].angle1[ia];
                                gamma = -1 * this->p_ucell_->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;

                                for(int m = 0; m < 2*l+1; m++)
                                {
                                    const int lm = l*l +m;
                                    if(index+2*l+1 > this->p_ucell_->natomwfc)
                                    {
                                        std::cout<<__FILE__<<__LINE__<<" "<<index<<" "<<this->p_ucell_->natomwfc<<std::endl;
                                        //ModuleBase::WARNING_QUIT("psi_initializer_atomic<T, Device>::cal_psig()","error: too many wfcs");
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
                                        (*(this->psig))(index, ig) = 
                                            this->template cast_to_T<T>(phase_factor(0.5*gamma, 0)*fup);
                                        (*(this->psig))(index, ig+this->pw_wfc_->npwk_max) = 
                                            this->template cast_to_T<T>(phase_factor(-0.5*gamma, 0)*fdw);
                                        //second rotation with angle gamma around(OZ)
                                        fup = phase_factor(0.5*(alpha + ModuleBase::PI),  1)*aux[ig];
                                        fdw = phase_factor(0.5*(alpha + ModuleBase::PI), -1)*aux[ig];
                                        (*(this->psig))(index+2*l+1, ig) = 
                                            this->template cast_to_T<T>(phase_factor(0.5*gamma, 0)*fup);
                                        (*(this->psig))(index+2*l+1, ig+this->pw_wfc_->npwk_max) = 
                                            this->template cast_to_T<T>(phase_factor(-0.5*gamma, 0)*fdw);
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
                                    //ModuleBase::WARNING_QUIT("psi_initializer_atomic<T, Device>::cal_psig()","error: too many wfcs");
                                }
                                for(int ig = 0;ig<npw;ig++)
                                {
                                     aux[ig] = sk[ig] * ylm(lm,ig) * ovlp_pswfcjlg[ig];
                                }
                                //rotate function
                                //first, rotation with angle alpha around(OX)
                                for(int ig = 0; ig<npw; ig++)
                                {
                                     fup = cos(0.5*alpha) * aux[ig];
                                     fdown = ModuleBase::IMAG_UNIT * sin(0.5* alpha) * aux[ig];
                                     //build the orthogonal wfc
                                     //first rotation with angle(alpha+ModuleBase::PI) around(OX)
                                     (*(this->psig))(index, ig) = 
                                        this->template cast_to_T<T>(
                                            (cos(0.5*gamman) + ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fup
                                                );
                                     (*(this->psig))(index, ig+ this->pw_wfc_->npwk_max) = 
                                        this->template cast_to_T<T>(
                                            (cos(0.5*gamman) - ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fdown
                                                );
                                     //second rotation with angle gamma around(OZ)
                                     fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                     fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                     (*(this->psig))(index+2*l+1, ig) = 
                                        this->template cast_to_T<T>(
                                            (cos(0.5*gamman) + ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fup
                                                );
                                     (*(this->psig))(index+2*l+1, ig+ this->pw_wfc_->npwk_max) = 
                                        this->template cast_to_T<T>(
                                            (cos(0.5*gamman) - ModuleBase::IMAG_UNIT * sin(0.5*gamman)) * fdown
                                                );
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
                                (*(this->psig))(index, ig) = 
                                    this->template cast_to_T<T>(
                                        lphase * sk [ig] * ylm(lm, ig) * ovlp_pswfcjlg[ig]
                                            );
                            }
                            index++;
                        }
                    }
                }
            }
			delete [] sk;
        }
    }

    delete[] ovlp_pswfcjlg;
    delete[] gk;
    delete[] aux;
    delete[] chiaux;
	/* complement the rest of bands if there are */
	if(this->nbands_complem() > 0)
	{
		this->random_t(this->psig->get_pointer(), index, this->psig->get_nbands(), ik);
	}
    ModuleBase::timer::tick("psi_initializer_atomic", "cal_psig");
    return this->psig;
}

template class psi_initializer_atomic<std::complex<double>, psi::DEVICE_CPU>;
template class psi_initializer_atomic<std::complex<float>, psi::DEVICE_CPU>;
// gamma point calculation
template class psi_initializer_atomic<double, psi::DEVICE_CPU>;
template class psi_initializer_atomic<float, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class psi_initializer_atomic<std::complex<double>, psi::DEVICE_GPU>;
template class psi_initializer_atomic<std::complex<float>, psi::DEVICE_GPU>;
// gamma point calculation
template class psi_initializer_atomic<double, psi::DEVICE_GPU>;
template class psi_initializer_atomic<float, psi::DEVICE_GPU>;
#endif