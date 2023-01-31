#include "wf_igk.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_psi/kernels/device.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/wf_op.h"

WF_igk::WF_igk()
{
}

WF_igk::~WF_igk()
{
	if (GlobalV::test_deconstructor)
	{
		std::cout << " ~WF_igk()" << std::endl;
	}
}


//--------------------------------------------------------
// Compute kinetic energy for each k-point
//--------------------------------------------------------

ModuleBase::Vector3<double> WF_igk::get_1qvec_cartesian(const int ik, const int ig) const
{
	ModuleBase::Vector3<double> qvec = GlobalC::wfcpw->getgpluskcar(ik,ig);

	/*
	if(igk(ik,ig)==0)
	{
		std::cout << " g add = " << GlobalC::sf.gcar << std::endl;
		std::cout << "\n igk = " << igk(ik,ig);
		std::cout << " k=" << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " <<
	GlobalC::kv.kvec_c[ik].z; std::cout << " g=" << GlobalC::sf.gcar[ this->igk(ik, ig) ].x << " " << GlobalC::sf.gcar[
	this->igk(ik, ig) ].y << " " << GlobalC::sf.gcar[ this->igk(ik, ig) ].z;
	}
	*/

	return qvec;
}

double *WF_igk::get_qvec_cartesian(const int &ik)
{
	double *qmod = new double[GlobalC::kv.ngk[ik]];
	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		// cartesian coordinate
		// modulus, in GlobalC::ucell.tpiba unit.
		const double q2 = GlobalC::wfcpw->getgk2(ik,ig);
		qmod[ig] = GlobalC::ucell.tpiba * sqrt(q2); // sqrt(q2) * GlobalC::ucell.tpiba;
	}
	return qmod;
}

std::complex<double> *WF_igk::get_sk(const int ik, const int it, const int ia, ModulePW::PW_Basis_K* wfc_basis) const
{
	ModuleBase::timer::tick("WF_igk", "get_sk");
	const double arg = (GlobalC::kv.kvec_c[ik] * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
	const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));
	std::complex<double> *sk = new std::complex<double>[GlobalC::kv.ngk[ik]];
	const int nx = wfc_basis->nx, ny = wfc_basis->ny, nz = wfc_basis->nz;
	for(int igl = 0; igl < GlobalC::kv.ngk[ik]; ++igl)
	{
		const int isz = wfc_basis->getigl2isz(ik,igl);
    	int iz = isz % nz;
    	const int is = isz / nz;
		const int ixy = wfc_basis->is2fftixy[is];
    	int ix = ixy / wfc_basis->fftny;
    	int iy = ixy % wfc_basis->fftny;
		if (ix < int(nx/2) + 1) ix += nx;
        if (iy < int(ny/2) + 1) iy += ny;
        if (iz < int(nz/2) + 1) iz += nz;
		const int iat = GlobalC::ucell.itia2iat(it, ia);
		sk[igl] = kphase * GlobalC::sf.eigts1(iat, ix) * GlobalC::sf.eigts2(iat, iy)
				 * GlobalC::sf.eigts3(iat, iz);
	}
	ModuleBase::timer::tick("WF_igk", "get_sk");
	return sk;
}

template<typename FPTYPE, typename Device>
void WF_igk::get_sk(Device * ctx, const int ik, ModulePW::PW_Basis_K* wfc_basis,  std::complex<FPTYPE> * sk) const
{
    ModuleBase::timer::tick("WF_igk", "get_sk");

    psi::DEVICE_CPU * cpu_ctx = {};
    psi::AbacusDevice_t device = psi::device::get_device_type<Device>(ctx);
    using cal_sk_op = hamilt::cal_sk_op<FPTYPE, Device>;
    using resmem_int_op = psi::memory::resize_memory_op<int, Device>;
    using delmem_int_op = psi::memory::delete_memory_op<int, Device>;
    using syncmem_int_op = psi::memory::synchronize_memory_op<int, Device, psi::DEVICE_CPU>;

    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;

    int iat = 0, _npw = GlobalC::kv.ngk[ik], eigts1_nc = GlobalC::sf.eigts1.nc, eigts2_nc = GlobalC::sf.eigts2.nc, eigts3_nc = GlobalC::sf.eigts3.nc;
    int * igl2isz = nullptr, * is2fftixy = nullptr, * atom_na = nullptr, * h_atom_na = new int[GlobalC::ucell.ntype];
    FPTYPE * atom_tau = nullptr, * h_atom_tau = new FPTYPE[GlobalC::ucell.nat * 3], * kvec = nullptr;
    std::complex<FPTYPE> * eigts1 = GlobalC::sf.get_eigts1_data<FPTYPE>(),
         * eigts2 = GlobalC::sf.get_eigts2_data<FPTYPE>(),
         * eigts3 = GlobalC::sf.get_eigts3_data<FPTYPE>();
    for (int it = 0; it < GlobalC::ucell.ntype; it++) {
        h_atom_na[it] = GlobalC::ucell.atoms[it].na;
        for (int ia = 0; ia < h_atom_na[it]; ia++) {
            auto * tau = reinterpret_cast<double *>(GlobalC::ucell.atoms[it].tau);
            h_atom_tau[iat * 3 + 0] = static_cast<FPTYPE>(tau[ia * 3 + 0]);
            h_atom_tau[iat * 3 + 1] = static_cast<FPTYPE>(tau[ia * 3 + 1]);
            h_atom_tau[iat * 3 + 2] = static_cast<FPTYPE>(tau[ia * 3 + 2]);
            iat++;
        }
    }
    if (device == psi::GpuDevice) {
        resmem_int_op()(ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_op()(ctx, cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);

        resmem_var_op()(ctx, atom_tau, GlobalC::ucell.nat * 3);
        resmem_var_op()(ctx, kvec, GlobalC::kv.kvec_c.size() * 3);
        syncmem_var_op()(ctx, cpu_ctx, atom_tau, h_atom_tau, GlobalC::ucell.nat * 3);
        syncmem_var_op()(ctx, cpu_ctx, kvec, reinterpret_cast<FPTYPE *>(GlobalC::kv.kvec_c.data()), GlobalC::kv.kvec_c.size() * 3);

        igl2isz = wfc_basis->d_igl2isz_k;
        is2fftixy = wfc_basis->d_is2fftixy;
    }
    else {
        atom_na = h_atom_na;
        atom_tau = h_atom_tau;
        igl2isz = wfc_basis->igl2isz_k;
        is2fftixy = wfc_basis->is2fftixy;
        kvec = reinterpret_cast<FPTYPE *>(GlobalC::kv.kvec_c.data());
    }

    cal_sk_op()(
        ctx,
        ik, GlobalC::ucell.ntype,
        wfc_basis->nx, wfc_basis->ny, wfc_basis->nz,
        _npw, wfc_basis->npwk_max,
        wfc_basis->fftny,
        eigts1_nc, eigts2_nc, eigts3_nc,
        atom_na, igl2isz, is2fftixy,
        ModuleBase::TWO_PI,
        kvec,
        atom_tau,
        eigts1, eigts2, eigts3,
        sk);
    if (device == psi::GpuDevice) {
        delmem_int_op()(ctx, atom_na);
        delmem_var_op()(ctx, kvec);
        delmem_var_op()(ctx, atom_tau);
    }
    delete [] h_atom_na;
    delete [] h_atom_tau;
    ModuleBase::timer::tick("WF_igk", "get_sk");
}

std::complex<double> *WF_igk::get_skq(int ik,
									  const int it,
									  const int ia,
									  ModuleBase::Vector3<double> q) // pengfei 2016-11-23
{
	std::complex<double> *skq = new std::complex<double>[GlobalC::kv.ngk[ik]];

	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		ModuleBase::Vector3<double> qkq = GlobalC::wfcpw->getgpluskcar(ik,ig) + q;
		double arg = (qkq * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
		skq[ig] = std::complex<double>(cos(arg), -sin(arg));
	}

	return skq;
}

template void WF_igk::get_sk<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int, ModulePW::PW_Basis_K*, std::complex<float>*) const;
template void WF_igk::get_sk<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int, ModulePW::PW_Basis_K*, std::complex<double>*) const;
#if defined(__CUDA) || defined(__ROCM)
template void WF_igk::get_sk<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int, ModulePW::PW_Basis_K*, std::complex<float>*) const;
template void WF_igk::get_sk<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int, ModulePW::PW_Basis_K*, std::complex<double>*) const;
#endif