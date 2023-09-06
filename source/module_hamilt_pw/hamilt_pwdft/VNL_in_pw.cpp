#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "VNL_in_pw.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/math_polyint.h"
#include "module_base/math_ylmreal.h"
#include "module_hamilt_pw/hamilt_pwdft/soc.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#include "module_psi/kernels/device.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/vnl_op.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            delmem_sd_op()(gpu_ctx, this->s_deeq);
            delmem_sd_op()(gpu_ctx, this->s_nhtol);
            delmem_sd_op()(gpu_ctx, this->s_nhtolm);
            delmem_sd_op()(gpu_ctx, this->s_indv);
            delmem_sd_op()(gpu_ctx, this->s_tab);
            delmem_cd_op()(gpu_ctx, this->c_deeq_nc);
            delmem_cd_op()(gpu_ctx, this->c_vkb);
        }
        else {
            delmem_zd_op()(gpu_ctx, this->z_deeq_nc);
        }
        delmem_dd_op()(gpu_ctx, this->d_deeq);
        delmem_zd_op()(gpu_ctx, this->z_vkb);
        delmem_dd_op()(gpu_ctx, this->d_tab);
        delmem_dd_op()(gpu_ctx, this->d_indv);
        delmem_dd_op()(gpu_ctx, this->d_nhtol);
        delmem_dd_op()(gpu_ctx, this->d_nhtolm);
    }
    else {
        if (GlobalV::precision_flag == "single") {
            delmem_sh_op()(cpu_ctx, this->s_deeq);
            delmem_sh_op()(cpu_ctx, this->s_nhtol);
            delmem_sh_op()(cpu_ctx, this->s_nhtolm);
            delmem_sh_op()(cpu_ctx, this->s_indv);
            delmem_sh_op()(cpu_ctx, this->s_tab);
            delmem_ch_op()(cpu_ctx, this->c_deeq_nc);
            delmem_ch_op()(cpu_ctx, this->c_vkb);
        }
        // There's no need to delete double precision pointers while in a CPU environment.
    }
}

//-----------------------------------
// setup lmaxkb, nhm, nkb, lmaxq 
// allocate vkb, GlobalV::NQX, tab, tab_at
//-----------------------------------
void pseudopot_cell_vnl::init(const int ntype,
                              Structure_Factor* psf_in,
                              const ModulePW::PW_Basis_K* wfc_basis,
                              const bool allocate_vkb)
{
	ModuleBase::TITLE("pseudopot_cell_vnl", "init");
	ModuleBase::timer::tick("ppcell_vnl", "init");

	GlobalV::ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << std::endl;

	int it = 0;
	this->wfcpw = wfc_basis;
    this->psf = psf_in;
    //----------------------------------------------------------
    // MEMBER VARIABLE :
    // NAME : lmaxkb(max angular momentum,(see pseudo_h))
    //----------------------------------------------------------
    this->lmaxkb = -1;
    for (it = 0;it < ntype; it++)
	{
		GlobalV::ofs_running << " " << GlobalC::ucell.atoms[it].label << " non-local projectors:" << std::endl;
		for (int ibeta = 0; ibeta < GlobalC::ucell.atoms[it].ncpp.nbeta; ibeta++) 
		{
			GlobalV::ofs_running << " projector " << ibeta+1 << " L=" << GlobalC::ucell.atoms[it].ncpp.lll[ibeta] <<  std::endl;
			this->lmaxkb = std::max( this->lmaxkb, GlobalC::ucell.atoms[it].ncpp.lll[ibeta]);
		}
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nhm(max number of different beta functions per atom)
//----------------------------------------------------------
	this->nhm = 0;
	for (it=0;it<ntype;it++)
	{	
		this->nhm = std::max(nhm, GlobalC::ucell.atoms[it].ncpp.nh);
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nkb(total number of beta functions, with struct.fact.)
//----------------------------------------------------------
	this->nkb = 0;
	for (it=0; it<ntype; it++)
	{
		this->nkb += GlobalC::ucell.atoms[it].ncpp.nh * GlobalC::ucell.atoms[it].na;
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"TOTAL NUMBER OF NONLOCAL PROJECTORS",nkb);

	if( this->nhm > 0 )
	{
		this->indv.create(ntype, this->nhm);
		this->nhtol.create(ntype, this->nhm);
		this->nhtolm.create(ntype, this->nhm);
		this->nhtoj.create(ntype, this->nhm);
		this->deeq.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
		this->deeq_nc.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm, this->nhm);
        if (GlobalV::device_flag == "gpu") {
            if (GlobalV::precision_flag == "single") {
                resmem_sd_op()(gpu_ctx, s_deeq, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
                resmem_sd_op()(gpu_ctx, s_nhtol, ntype * this->nhm);
                resmem_sd_op()(gpu_ctx, s_nhtolm, ntype * this->nhm);
                resmem_sd_op()(gpu_ctx, s_indv, ntype * this->nhm);
                resmem_cd_op()(gpu_ctx, c_deeq_nc, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            }
            else {
                resmem_zd_op()(gpu_ctx, z_deeq_nc, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            }
            resmem_dd_op()(gpu_ctx, d_deeq, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            resmem_dd_op()(gpu_ctx, d_indv, ntype * this->nhm);
            resmem_dd_op()(gpu_ctx, d_nhtol, ntype * this->nhm);
            resmem_dd_op()(gpu_ctx, d_nhtolm, ntype * this->nhm);
        }
        else {
            if (GlobalV::precision_flag == "single") {
                resmem_sh_op()(cpu_ctx, s_deeq, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm, "VNL::s_deeq");
                resmem_sh_op()(cpu_ctx, s_nhtol, ntype * this->nhm, "VNL::s_nhtol");
                resmem_sh_op()(cpu_ctx, s_nhtolm, ntype * this->nhm, "VNL::s_nhtolm");
                resmem_sh_op()(cpu_ctx, s_indv, ntype * this->nhm, "VNL::s_indv");
                resmem_ch_op()(cpu_ctx, c_deeq_nc, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm, "VNL::c_deeq_nc");
            }
            else {
                this->z_deeq_nc = this->deeq_nc.ptr;
            }
            this->d_deeq = this->deeq.ptr;
            this->d_indv = this->indv.c;
            this->d_nhtol = this->nhtol.c;
            this->d_nhtolm = this->nhtolm.c;
            // There's no need to delete double precision pointers while in a CPU environment.
        }
		this->dvan.create(ntype, this->nhm, this->nhm);
		this->dvan_so.create(GlobalV::NSPIN, ntype, this->nhm, this->nhm);
		this->becsum.create(GlobalV::NSPIN, GlobalC::ucell.nat, this->nhm * (this->nhm + 1) / 2);
	}
	else
	{
		GlobalV::ofs_running << "\n nhm = 0, not allocate some matrix.";
	}

	// nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
	// dq+4)*cell_factor;
	this->lmaxq = 2 * this->lmaxkb + 1;
	int npwx = this->wfcpw->npwk_max;
	if (nkb > 0 && allocate_vkb )
	{
		vkb.create(nkb, npwx);
		ModuleBase::Memory::record("VNL::vkb", nkb * npwx * sizeof(double));
	}

	//this->nqx = 10000;		// calculted in allocate_nlpot.f90
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ); //LiuXh modify 20180515
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) + 1000; //LiuXh add 20180515
	//GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) * 10; //LiuXh add 20180515
	GlobalV::NQX = this->calculate_nqx(INPUT.ecutwfc,GlobalV::DQ) * cell_factor; //LiuXh add 20180619
	// nqx = (sqrt(ecutwfc)/dq+4)*cell_factor;

	
	// mohan update 2021-02-22
	const int nbrx = 10;
	const int nbrx_nc = 20;
	//  max number of beta functions
	if(GlobalV::NSPIN!=4) 
	{
		this->tab.create(ntype, nbrx, GlobalV::NQX);
		ModuleBase::Memory::record("VNL::tab", ntype * nbrx * GlobalV::NQX * sizeof(double));
	}
	else 
	{
		this->tab.create(ntype, nbrx_nc, GlobalV::NQX);
		ModuleBase::Memory::record("VNL::tab", ntype * nbrx_nc * GlobalV::NQX * sizeof(double));
	}

	
	// mohan update 2021-02-22
	int nchix = 10;
	int nchix_nc = 20;
	// nchix : max number of atomic wavefunctions per atom
	if(GlobalV::NSPIN!=4) 
	{
		this->tab_at.create(ntype, nchix, GlobalV::NQX);
		ModuleBase::Memory::record("VNL::tab_at", ntype * nchix * GlobalV::NQX * sizeof(double));
	}
	else 
	{
		this->tab_at.create(ntype, nchix_nc, GlobalV::NQX);
		ModuleBase::Memory::record("VNL::tab_at", ntype * nchix_nc * GlobalV::NQX * sizeof(double));
	}
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            resmem_sd_op()(gpu_ctx, s_tab, this->tab.getSize());
            resmem_cd_op()(gpu_ctx, c_vkb, nkb * npwx);
        }
        resmem_zd_op()(gpu_ctx, z_vkb, nkb * npwx);
        resmem_dd_op()(gpu_ctx, d_tab, this->tab.getSize());
    }
    else {
        if (GlobalV::precision_flag == "single") {
            resmem_sh_op()(cpu_ctx, s_tab, this->tab.getSize());
            resmem_ch_op()(cpu_ctx, c_vkb, nkb * npwx);
        }
        this->z_vkb = this->vkb.c;
        this->d_tab = this->tab.ptr;
        // There's no need to delete double precision pointers while in a CPU environment.
    }

	ModuleBase::timer::tick("ppcell_vnl","init");
	return;
}



//----------------------------------------------------------
// Calculates beta functions (Kleinman-Bylander projectors),
// with structure factor, for all atoms, in reciprocal space
//----------------------------------------------------------
void pseudopot_cell_vnl::getvnl(const int &ik, ModuleBase::ComplexMatrix& vkb_in)const
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl");
	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	if(lmaxkb < 0) 
	{
		return;
	}

	const int npw = this->wfcpw->npwk[ik];

	// When the internal memory is large enough, it is better to make vkb1 be the number of pseudopot_cell_vnl.
    // We only need to initialize it once as long as the cell is unchanged.
	ModuleBase::matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix ylm(x1, npw);
	ModuleBase::Memory::record("VNL::ylm", x1 * npw * sizeof(double));
	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
	for (int ig = 0;ig < npw;ig++) 
	{
		gk[ig] = this->wfcpw->getgpluskcar(ik, ig);
	}

	ModuleBase::YlmReal::Ylm_Real(cpu_ctx, x1, npw, reinterpret_cast<double *>(gk), ylm.c);

    using Device = psi::DEVICE_CPU;
    Device * ctx = {};
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<double>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<double>, Device>;
    std::complex<double> * sk = nullptr;
    resmem_complex_op()(ctx, sk, GlobalC::ucell.nat * npw, "VNL::sk");
    this->psf->get_sk(ctx, ik, this->wfcpw, sk);

    int jkb = 0, iat = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (int nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (int ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;

				vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
						this->tab, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );
			}

			// add spherical harmonic part
			for (int ih = 0;ih < nh;ih++)
			{
				if (nb == this->indv(it, ih))
				{
					const int lm = static_cast<int>( nhtolm(it, ih) );
					for (int ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++) 
		{
			for (int ih = 0;ih < nh;ih++)
			{
				std::complex<double> pref = pow( ModuleBase::NEG_IMAG_UNIT, nhtol(it, ih));	//?
				std::complex<double>* pvkb = &vkb_in(jkb, 0);
				for (int ig = 0;ig < npw;ig++)
				{
					pvkb[ig] = vkb1(ih, ig) * sk [iat * npw + ig] * pref;
				}
				++jkb;
			} // end ih
            iat ++;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;
    delmem_complex_op()(ctx, sk);
	ModuleBase::timer::tick("pp_cell_vnl","getvnl");

	return;
} // end subroutine getvnl

template <typename FPTYPE, typename Device>
void pseudopot_cell_vnl::getvnl(Device * ctx, const int &ik, std::complex<FPTYPE>* vkb_in)const
{
    if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl");
    ModuleBase::timer::tick("pp_cell_vnl","getvnl");

    using cal_vnl_op = hamilt::cal_vnl_op<FPTYPE, Device>;
    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using castmem_var_h2d_op = psi::memory::cast_memory_op<FPTYPE, double, Device, psi::DEVICE_CPU>;
    using castmem_var_h2h_op = psi::memory::cast_memory_op<FPTYPE, double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;

    if(lmaxkb < 0)
    {
        return;
    }

    const int x1 = (lmaxkb + 1) * (lmaxkb + 1);
    const int npw = this->wfcpw->npwk[ik];

    int * atom_nh = nullptr, * atom_na = nullptr, * atom_nb = nullptr, * h_atom_nh = new int[GlobalC::ucell.ntype], * h_atom_na = new int[GlobalC::ucell.ntype], * h_atom_nb = new int[GlobalC::ucell.ntype];
    for (int it = 0; it < GlobalC::ucell.ntype; it++) {
        h_atom_nb[it] = GlobalC::ucell.atoms[it].ncpp.nbeta;
        h_atom_nh[it] = GlobalC::ucell.atoms[it].ncpp.nh;
        h_atom_na[it] = GlobalC::ucell.atoms[it].na;
    }
    // When the internal memory is large enough, it is better to make vkb1 be the number of pseudopot_cell_vnl.
    // We only need to initialize it once as long as the cell is unchanged.
    FPTYPE * vkb1 = nullptr, * gk = nullptr, * ylm = nullptr,
           * _tab = this->get_tab_data<FPTYPE>(),
           * _indv = this->get_indv_data<FPTYPE>(),
           * _nhtol = this->get_nhtol_data<FPTYPE>(),
           * _nhtolm = this->get_nhtolm_data<FPTYPE>();
    resmem_var_op()(ctx, ylm, x1 * npw, "VNL::ylm");
    resmem_var_op()(ctx, vkb1, nhm * npw, "VNL::vkb1");

    ModuleBase::Vector3<double> *_gk = new ModuleBase::Vector3<double>[npw];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
    for (int ig = 0;ig < npw; ig++)
    {
        _gk[ig] = this->wfcpw->getgpluskcar(ik,ig);
    }
    if (GlobalV::device_flag == "gpu") {
        resmem_int_op()(ctx, atom_nh, GlobalC::ucell.ntype);
        resmem_int_op()(ctx, atom_nb, GlobalC::ucell.ntype);
        resmem_int_op()(ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_nh, h_atom_nh, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_nb, h_atom_nb, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);

        resmem_var_op()(ctx, gk, npw * 3);
        castmem_var_h2d_op()(ctx, cpu_ctx, gk, reinterpret_cast<double *>(_gk), npw * 3);
    }
    else {
        atom_nh = h_atom_nh;
        atom_nb = h_atom_nb;
        atom_na = h_atom_na;
        if (GlobalV::precision_flag == "single") {
            resmem_var_op()(ctx, gk, npw * 3);
            castmem_var_h2h_op()(cpu_ctx, cpu_ctx, gk, reinterpret_cast<double *>(_gk), npw * 3);
        }
        else {
            gk = reinterpret_cast<FPTYPE *>(_gk);
        }
    }

    ModuleBase::YlmReal::Ylm_Real(ctx, x1, npw, gk, ylm);

    std::complex<FPTYPE> * sk = nullptr;
    resmem_complex_op()(ctx, sk, GlobalC::ucell.nat * npw);
    this->psf->get_sk(ctx, ik, this->wfcpw, sk);

#ifdef USE_PAW
	if(GlobalV::use_paw)
	{
		double* kpt;
		kpt = new double[3];
		kpt[0] = this->wfcpw->kvec_c[ik].x;
		kpt[1] = this->wfcpw->kvec_c[ik].y;
		kpt[2] = this->wfcpw->kvec_c[ik].z;

		double ** kpg;
		kpg = new double*[npw];
		for(int ipw=0;ipw<npw;ipw++)
		{
			kpg[ipw] = new double[3];
			kpg[ipw][0] = _gk[ipw].x;
			kpg[ipw][1] = _gk[ipw].y;
			kpg[ipw][2] = _gk[ipw].z;
		}

		GlobalC::paw_cell.set_paw_k(npw,kpt,
			this->wfcpw->get_ig2ix(ik).data(),
			this->wfcpw->get_ig2iy(ik).data(),
			this->wfcpw->get_ig2iz(ik).data(),
			(const double **) kpg,GlobalC::ucell.tpiba);

		delete[] kpt;
		for(int ipw = 0; ipw < npw; ipw++)
		{
			delete[] kpg[ipw];
		}
		delete[] kpg;

		GlobalC::paw_cell.get_vkb();
	}
#endif

    cal_vnl_op()(
        ctx,
        GlobalC::ucell.ntype,  npw, this->wfcpw->npwk_max, this->nhm, GlobalV::NQX,
        this->tab.getBound2(), this->tab.getBound3(),
        atom_na, atom_nb, atom_nh,
        static_cast<FPTYPE>(GlobalV::DQ), static_cast<FPTYPE>(GlobalC::ucell.tpiba), static_cast<std::complex<FPTYPE>>(ModuleBase::NEG_IMAG_UNIT),
        gk, ylm, _indv, _nhtol, _nhtolm, _tab, vkb1, sk,
        vkb_in);

    delete [] _gk;
    delete [] h_atom_nh;
    delete [] h_atom_na;
    delete [] h_atom_nb;
    delmem_var_op()(ctx, ylm);
    delmem_var_op()(ctx, vkb1);
    delmem_complex_op()(ctx, sk);
    if (psi::device::get_device_type<Device>(ctx) == psi::GpuDevice) {
        delmem_var_op()(ctx, gk);
        delmem_int_op()(ctx, atom_nh);
        delmem_int_op()(ctx, atom_nb);
        delmem_int_op()(ctx, atom_na);
    }
    ModuleBase::timer::tick("pp_cell_vnl","getvnl");
} // end subroutine getvnl

void pseudopot_cell_vnl::init_vnl(UnitCell &cell)
{
	ModuleBase::TITLE("pseudopot_cell_vnl","init_vnl");
	ModuleBase::timer::tick("ppcell_vnl","init_vnl");

	//from init_us_1
	//   a) For each non vanderbilt pseudopotential it computes the D and
	//      the betar in the same form of the Vanderbilt pseudopotential.
	//   b) It computes the indices indv which establish the correspondence
	//      nh <-> beta in the atom
	//   c) It computes the indices nhtol which establish the correspondence
	//      nh <-> angular momentum of the beta function
	//   d) It computes the indices nhtolm which establish the correspondence
	//      nh <-> combined (l,m) index for the beta function.

	// For each pseudopotential we initialize the indices nhtol, nhtolm,
	// nhtoj, indv, and if the pseudopotential is of KB type we initialize
	// the atomic D terms
	this->dvan.zero_out();
	this->dvan_so.zero_out();//added by zhengdy-soc

	for(int it=0;it<cell.ntype;it++)
	{
		int BetaIndex=0;
		const int Nprojectors = cell.atoms[it].ncpp.nh;
		for (int ib=0; ib<cell.atoms[it].ncpp.nbeta; ib++)
		{
			const int l = cell.atoms[it].ncpp.lll [ib];
			const double j = cell.atoms[it].ncpp.jjj [ib];
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
				this->nhtoj(it,BetaIndex) = j;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}

		//    From now on the only difference between KB and US pseudopotentials
		//    is in the presence of the q and Q functions.
		//    Here we initialize the D of the solid
		if(cell.atoms[it].ncpp.has_so )
		{
			Soc soc;
			soc.rot_ylm(this->lmaxkb);
			soc.fcoef.create(cell.ntype, this->nhm, this->nhm);
			for(int ip=0; ip<Nprojectors; ip++)
			{
				const int l1 = this->nhtol (it, ip);
				const double j1 = this->nhtoj (it, ip);
				const int m1 = this->nhtolm (it, ip) - l1* l1;
				//const int v1 = static_cast<int>( indv(it, ip ) );
				for(int ip2=0;ip2<Nprojectors; ip2++)
				{
					const int l2 = this->nhtol (it, ip2);
					const double j2 = this->nhtoj (it, ip2);
					const int m2 = this->nhtolm (it, ip2) - l2* l2;
					//const int v2 = static_cast<int>( indv(it, ip2 ) );
					if(l1 == l2 && fabs(j1-j2)<1e-7)
					{
						for(int is1=0;is1<2;is1++)
						{
							for(int is2=0;is2<2;is2++)
							{
								soc.set_fcoef(l1, l2,
										is1, is2,
										m1, m2,
										j1, j2,
										it, ip, ip2);
							}
						}
					}
				}
			}
//
//   and calculate the bare coefficients
//
			for(int ip = 0;ip<Nprojectors; ++ip)
			{
				const int ir = static_cast<int>( indv(it, ip ) );
				for(int ip2=0; ip2<Nprojectors; ++ip2)
				{
					const int is = static_cast<int>( indv(it, ip2) );
					int ijs =0;
					for(int is1=0;is1<2;++is1)
					{
						for(int is2=0;is2<2;++is2)
						{
							this->dvan_so(ijs,it,ip,ip2) = cell.atoms[it].ncpp.dion(ir, is) * soc.fcoef(it,is1,is2,ip,ip2);
							++ijs;
							if(ir != is) soc.fcoef(it,is1,is2,ip,ip2) = std::complex<double>(0.0,0.0);
						}
					}
				}
			}
		}
		else
		for (int ip=0; ip<Nprojectors; ip++)
		{
			for (int ip2=0; ip2<Nprojectors; ip2++)
			{
				if ( this->nhtol (it, ip) == nhtol (it, ip2) &&
				     this->nhtolm(it, ip) == nhtolm(it, ip2) )
				{
					const int ir = static_cast<int>( indv(it, ip ) );
					const int is = static_cast<int>( indv(it, ip2) );
					if(GlobalV::LSPINORB)
					{
						this->dvan_so(0,it,ip,ip2) = cell.atoms[it].ncpp.dion(ir, is);
						this->dvan_so(3,it,ip,ip2) = cell.atoms[it].ncpp.dion(ir, is);
					}
					else
					{
						this->dvan(it, ip, ip2) = cell.atoms[it].ncpp.dion(ir, is);
					}
				}
			} 
		} 
	} 

	// h) It fills the interpolation table for the beta functions
	/**********************************************************
	// He Lixin: this block is used for non-local potential
	// fill the interpolation table tab
	************************************************************/

	const double pref = ModuleBase::FOUR_PI / sqrt(cell.omega);
	this->tab.zero_out();
	GlobalV::ofs_running<<"\n Init Non-Local PseudoPotential table : ";
	for (int it = 0;it < cell.ntype;it++)  
	{
		const int nbeta = cell.atoms[it].ncpp.nbeta;
		int kkbeta = cell.atoms[it].ncpp.kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			const int l = cell.atoms[it].ncpp.lll[ib];
			for (int iq=0; iq<GlobalV::NQX; iq++)  
			{
				const double q = iq * GlobalV::DQ;
				ModuleBase::Sphbes::Spherical_Bessel(kkbeta, cell.atoms[it].ncpp.r, q, l, jl);

				for (int ir = 0;ir < kkbeta;ir++)
				{
					aux[ir] = cell.atoms[it].ncpp.betar(ib, ir) *
					          jl[ir] * cell.atoms[it].ncpp.r[ir];
				} 
				double vqint;
				ModuleBase::Integral::Simpson_Integral(kkbeta, aux, cell.atoms[it].ncpp.rab, vqint);
				this->tab(it, ib, iq) = vqint * pref;
			} 
		} 
		delete[] aux;
		delete[] jl;
	}
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_tab, this->tab.ptr, this->tab.getSize());
        }
        // Even when the single precision flag is enabled,
        // these variables are utilized in the Force/Stress calculation as well.
        // modified by denghuilu at 2023-05-15
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_indv, this->indv.c, this->indv.nr * this->indv.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_tab, this->tab.ptr, this->tab.getSize());
    }
    else {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_indv, this->indv.c, this->indv.nr * this->indv.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_nhtol, this->nhtol.c, this->nhtol.nr * this->nhtol.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_nhtolm, this->nhtolm.c, this->nhtolm.nr * this->nhtolm.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_tab, this->tab.ptr, this->tab.getSize());
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
	ModuleBase::timer::tick("ppcell_vnl","init_vnl");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done." << std::endl;
	return;
}

#ifdef __LCAO
std::complex<double> pseudopot_cell_vnl::Cal_C(int alpha, int lu, int mu, int L, int M)   // pengfei Li  2018-3-23
{
	std::complex<double> cf;
	if(alpha == 0)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,1,L,M);
	}
	else if(alpha == 1)
	{
		cf = -sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,2,L,M);
	}
	else if(alpha == 2)
	{
		cf = sqrt(4*ModuleBase::PI/3)*CG(lu,mu,1,0,L,M);
	}
	else
	{
		ModuleBase::WARNING_QUIT("pseudopot_cell_vnl_alpha", "alpha must be 0~2");
	}
	
	return cf;
}

double pseudopot_cell_vnl::CG(int l1, int m1, int l2, int m2, int L, int M)      // pengfei Li 2018-3-23
{
	int dim = L*L+M;
	int dim1 = l1*l1+m1;
	int dim2 = l2*l2+m2;
	
	//double A = MGT.Gaunt_Coefficients(dim1, dim2, dim);
	
	return MGT.Gaunt_Coefficients(dim1, dim2, dim);
}

// void pseudopot_cell_vnl::getvnl_alpha(const int &ik)           // pengfei Li  2018-3-23
// {
// 	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","getvnl_alpha");
// 	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");

// 	if(lmaxkb < 0) 
// 	{
// 		return;
// 	}
	
// 	const int npw = this->wfcpw->npwk[ik];
// 	int ig, ia, nb, ih, lu, mu;

// 	double *vq = new double[npw];
// 	const int x1= (lmaxkb + 2)*(lmaxkb + 2);

// 	ModuleBase::matrix ylm(x1, npw);
// 	ModuleBase::Vector3<double> *gk = new ModuleBase::Vector3<double>[npw];
// 	for (ig = 0;ig < npw;ig++) 
// 	{
// 		gk[ig] = this->wfcpw->getgpluskcar(ik,ig);
// 	}

// 	vkb1_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb1_alpha[i] = new std::complex<double>*[nhm];
// 		for(int j=0; j<nhm; j++)
// 		{
// 			vkb1_alpha[i][j] = new std::complex<double>[npw];
// 		}
// 	}	
	
// 	vkb_alpha = new std::complex<double>**[3];
// 	for(int i=0; i<3; i++)
// 	{
// 		vkb_alpha[i] = new std::complex<double>*[nkb];
// 		for(int j=0; j<nkb; j++)
// 		{
// 			vkb_alpha[i][j] = new std::complex<double>[this->wfcpw->npwk_max];
// 		}
// 	}
	
// 	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

// 	MGT.init_Gaunt_CH( lmaxkb + 2 );
// 	MGT.init_Gaunt( lmaxkb + 2 );
	
// 	int jkb = 0;
// 	for(int it = 0;it < GlobalC::ucell.ntype;it++)
// 	{
// 		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
// 		// calculate beta in G-space using an interpolation table
// 		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
// 		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

// 		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

// 		for(int i=0; i<3; i++)
// 			for(int j=0; j<nhm; j++)
// 			{
// 				ModuleBase::GlobalFunc::ZEROS(vkb1_alpha[i][j], npw);
// 			}
			
// 		for (ih = 0;ih < nh; ih++)
// 		{
// 			lu = static_cast<int>( nhtol(it, ih));
// 			mu = static_cast<int>( nhtolm(it, ih)) - lu * lu;
// 			nb = static_cast<int>( indv(it, ih));
			
// 			for (int L= std::abs(lu - 1); L<= (lu + 1); L++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					const double gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
// 					vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
// 							this->tab_alpha, it, nb, L, GlobalV::NQX, GlobalV::DQ, gnorm);
					
// 					for (int M=0; M<2*L+1; M++)
// 					{
// 						int lm = L*L + M;
// 						for (int alpha=0; alpha<3; alpha++)
// 						{
// 							std::complex<double> c = Cal_C(alpha,lu, mu, L, M);
// 							/*if(alpha == 0)
// 							{
// 								std::cout<<"lu= "<<lu<<"  mu= "<<mu<<"  L= "<<L<<"  M= "<<M<<" alpha = "<<alpha<<"  "<<c<<std::endl;
// 							}*/
// 							vkb1_alpha[alpha][ih][ig] += c * vq[ig] * ylm(lm, ig) * pow( ModuleBase::NEG_IMAG_UNIT, L);
// 						}	
// 					}
// 				}
// 			}
// 		} // end nbeta

// 		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
// 		{
// 			std::complex<double> *sk = this->psf->get_sk(ik, it, ia,this->wfcpw);
// 			for (ih = 0;ih < nh;ih++)
// 			{
// 				for (ig = 0;ig < npw;ig++)
// 				{
// 					for(int alpha=0; alpha<3; alpha++)
// 					{
// 						vkb_alpha[alpha][jkb][ig] = vkb1_alpha[alpha][ih][ig] * sk [ig];
// 					}
// 				}
// 				++jkb;
// 			} // end ih
// 			delete [] sk;
// 		} // end ia
// 	} // enddo

// 	delete [] gk;
// 	delete [] vq;
// 	ModuleBase::timer::tick("pp_cell_vnl","getvnl_alpha");
// 	return;
// } 
#endif

void pseudopot_cell_vnl::init_vnl_alpha(void)          // pengfei Li 2018-3-23
{
	if(GlobalV::test_pp) ModuleBase::TITLE("pseudopot_cell_vnl","init_vnl_alpha");
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");

	for(int it=0;it<GlobalC::ucell.ntype;it++)
	{
		int BetaIndex=0;
		//const int Nprojectors = GlobalC::ucell.atoms[it].nh;
		for (int ib=0; ib<GlobalC::ucell.atoms[it].ncpp.nbeta; ib++)
		{
			const int l = GlobalC::ucell.atoms[it].ncpp.lll [ib];
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}
	} 


	// max number of beta functions
	const int nbrx = 10;

	const double pref = ModuleBase::FOUR_PI / sqrt(GlobalC::ucell.omega);
	this->tab_alpha.create(GlobalC::ucell.ntype, nbrx, lmaxkb+2, GlobalV::NQX);
	this->tab_alpha.zero_out();
	GlobalV::ofs_running<<"\n Init Non-Local PseudoPotential table( including L index) : ";
	for (int it = 0;it < GlobalC::ucell.ntype;it++)  
	{
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		int kkbeta = GlobalC::ucell.atoms[it].ncpp.kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			for (int L = 0; L <= lmaxkb+1; L++)
			{
				for (int iq = 0; iq < GlobalV::NQX; iq++)
				{
					const double q = iq * GlobalV::DQ;
					ModuleBase::Sphbes::Spherical_Bessel(kkbeta, GlobalC::ucell.atoms[it].ncpp.r, q, L, jl);
					
					for (int ir = 0;ir < kkbeta;ir++)
					{
						aux[ir] = GlobalC::ucell.atoms[it].ncpp.betar(ib, ir) * jl[ir] * 
								  GlobalC::ucell.atoms[it].ncpp.r[ir] * GlobalC::ucell.atoms[it].ncpp.r[ir];
					}
					double vqint;
					ModuleBase::Integral::Simpson_Integral(kkbeta, aux, GlobalC::ucell.atoms[it].ncpp.rab, vqint);
					this->tab_alpha(it, ib, L, iq) = vqint * pref;
				}
			}
		} 
		delete[] aux;
		delete[] jl;
	}
	ModuleBase::timer::tick("ppcell_vnl","init_vnl_alpha");
	GlobalV::ofs_running << "\n Init Non-Local-Pseudopotential done(including L)." << std::endl;
	return;
}



void pseudopot_cell_vnl::print_vnl(std::ofstream &ofs)
{
	output::printr3_d(ofs, " tab : ", tab);
}



int pseudopot_cell_vnl::calculate_nqx(const double &ecutwfc,const double &dq)
{
	int points_of_table = static_cast<int>( sqrt(ecutwfc)/dq + 4 ) ;//* cell_factor;
//	std::cout<<"\n nqx = "<< points_of_table;
//----------------------------------------------------------
// EXPLAIN : Plus 1 because the formula is transfered from 
// fortran code.
//----------------------------------------------------------
	return points_of_table + 1; 
}

// ----------------------------------------------------------------------
void pseudopot_cell_vnl::cal_effective_D(void)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "cal_effective_D");

    /*
	recalculate effective coefficient matrix for non-local pseudo-potential
	1. assign to each atom from element;
	2. extend to each spin when nspin larger than 1
	3. rotate to effective matrix when spin-orbital coupling is used
	*/

    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int nht = GlobalC::ucell.atoms[it].ncpp.nh;
        // nht: number of beta functions per atom type
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            for (int ih = 0; ih < nht; ih++)
            {
                for (int jh = ih; jh < nht; jh++)
                {
                    if (GlobalV::LSPINORB)
                    {
                        this->deeq_nc(is, iat, ih, jh) = this->dvan_so(is, it, ih, jh);
                        this->deeq_nc(is, iat, jh, ih) = this->dvan_so(is, it, jh, ih);
                    }
                    else if (GlobalV::NSPIN == 4)
                    {
                        if (is == 0)
                        {
                            this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                            this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                        }
                        else if (is == 1)
                        {
                            this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                            this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                        }
                        else if (is == 2)
                        {
                            this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                            this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                        }
                        else if (is == 3)
                        {
                            this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                            this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                        }
                    }
                    else
                    {
                        this->deeq(is, iat, ih, jh) = this->dvan(it, ih, jh);
                        this->deeq(is, iat, jh, ih) = this->dvan(it, ih, jh);
                        // in most of pseudopotential files, number of projections of one orbital is only one, 
                        // which lead to diagonal matrix of dion
                        // when number larger than 1, non-diagonal dion should be calculated.
                        if(ih != jh && std::fabs(this->deeq(is, iat, ih, jh))>0.0)
                        {
                            this->multi_proj = true;
                        }
                    }
                }
            }
        }
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, this->s_deeq, this->deeq.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, this->c_deeq_nc, this->deeq_nc.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
        }
        else {
            syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, this->z_deeq_nc, this->deeq_nc.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
        }
        syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, this->d_deeq, this->deeq.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
    }
    else {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, this->s_deeq, this->deeq.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
            castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, this->c_deeq_nc, this->deeq_nc.ptr, GlobalV::NSPIN * GlobalC::ucell.nat * this->nhm * this->nhm);
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
}

template <>
float * pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->s_nhtol;
}
template <>
double * pseudopot_cell_vnl::get_nhtol_data() const
{
    return this->d_nhtol;
}

template <>
float * pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->s_nhtolm;
}
template <>
double * pseudopot_cell_vnl::get_nhtolm_data() const
{
    return this->d_nhtolm;
}

template <>
float * pseudopot_cell_vnl::get_indv_data() const
{
    return this->s_indv;
}
template <>
double * pseudopot_cell_vnl::get_indv_data() const
{
    return this->d_indv;
}

template <>
float * pseudopot_cell_vnl::get_tab_data() const
{
    return this->s_tab;
}
template <>
double * pseudopot_cell_vnl::get_tab_data() const
{
    return this->d_tab;
}

template <>
float * pseudopot_cell_vnl::get_deeq_data() const
{
    return this->s_deeq;
}
template <>
double * pseudopot_cell_vnl::get_deeq_data() const
{
    return this->d_deeq;
}

template <>
std::complex<float> * pseudopot_cell_vnl::get_vkb_data() const
{
    return this->c_vkb;
}
template <>
std::complex<double> * pseudopot_cell_vnl::get_vkb_data() const
{
    return this->z_vkb;
}

template <>
std::complex<float> * pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->c_deeq_nc;
}
template <>
std::complex<double> * pseudopot_cell_vnl::get_deeq_nc_data() const
{
    return this->z_deeq_nc;
}

template void pseudopot_cell_vnl::getvnl<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int const&, std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int const&, std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
template void pseudopot_cell_vnl::getvnl<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int const&, std::complex<float>*) const;
template void pseudopot_cell_vnl::getvnl<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int const&, std::complex<double>*) const;
#endif
