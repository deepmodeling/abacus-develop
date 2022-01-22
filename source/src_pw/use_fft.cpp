#include "use_fft.h"
#include "global.h"
#include "../module_base/memory.h"

std::vector<std::complex<double>> Use_FFT::porter;

Use_FFT::Use_FFT()
{
}

Use_FFT::~Use_FFT()
{
#ifdef __CUDA
    cudaFree(d_porter);
    cufftDestroy(fft_handle);
#endif
}

std::complex<double>* Use_FFT::get_porter(const int &begin, const int &end)
{
    assert(begin >= 0);
    assert(end >= begin);
    assert(porter.size() >= end);
	return &(porter[begin]);
}

void Use_FFT::allocate(void)
{
    ModuleBase::TITLE("Use_FFT","allocate");

	this->porter.resize(GlobalC::pw.nrxx * GlobalV::NPOL);
#ifdef __CUDA
    cufftPlan3d(&fft_handle, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, CUFFT_Z2Z);
    cudaMalloc((void**)&d_porter, GlobalC::pw.nrxx * sizeof(double2));
#endif
    ModuleBase::Memory::record("Use_FFT","porter",GlobalC::pw.nrxx,"complexmatrix");
    return;
}

void Use_FFT::RoundTrip(
    const std::complex<double> *psi,
    const double *vr,
    const int *fft_index,
    const int &max_g,
    const int &max_r,
    std::complex<double> *psic)
{
	// (1) set value
	for (int ig=0; ig< max_g; ig++)
    {
        psic[ fft_index[ig]  ] = psi[ig];
    }

	// (2) fft to real space and doing things.
    GlobalC::pw.FFT_wfc.FFT3D( psic, 1);
    for (int ir=0; ir< max_r; ir++)
    {
        psic[ir] *= vr[ir];
    }

	// (3) fft back to G space.
    GlobalC::pw.FFT_wfc.FFT3D( psic, -1);
    return;
}

void Use_FFT::ToRealSpace(const int &is, const ModuleBase::ComplexMatrix &vg, double *vr)
{
	// (1) set value
    ModuleBase::GlobalFunc::ZEROS( porter.data(), GlobalC::pw.nrxx );
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter[ GlobalC::pw.ig2fftc[ig] ] = vg(is, ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_chg.FFT3D(porter.data(), 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToRealSpace_psi(const int &ik, const std::complex<double> *psig, std::complex<double> *psir)
{
	ModuleBase::GlobalFunc::ZEROS(psir, GlobalC::pw.nrxx);
	for(int ig=0; ig<GlobalC::wf.npw; ig++)
	{
		psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = psig[ig];
	}
	// 1: to real space.
	GlobalC::pw.FFT_wfc.FFT3D(psir, 1);

	return;
}


void Use_FFT::ToRealSpace_psi(const int &ik, const int &ib, const ModuleBase::ComplexMatrix &evc, std::complex<double> *psir)
{
	// (1) set value
    ModuleBase::GlobalFunc::ZEROS( psir, GlobalC::pw.nrxx );
    for (int ig=0; ig<GlobalC::wf.npw; ig++)
    {
        psir[ GlobalC::pw.ig2fftc[ GlobalC::wf.igk(ik,ig) ] ] = evc(ib, ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_wfc.FFT3D(psir, 1);
    return;
}



void Use_FFT::ToRealSpace(const int &is, const ModuleBase::ComplexMatrix &vg, ModuleBase::matrix &vr)
{
	// (1) set value
    ModuleBase::GlobalFunc::ZEROS( porter.data(), GlobalC::pw.nrxx);
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter [GlobalC::pw.ig2fftc[ig]] = vg(is,ig);
    }

	// (2) fft and get value
    GlobalC::pw.FFT_chg.FFT3D(porter.data(), 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr(is,ir) = porter[ir].real();
    }
    return;
}


// Fourer transform of vg,
// then put vg into vr.
void Use_FFT::ToRealSpace(const std::complex<double> *vg, double *vr)
{
    ModuleBase::GlobalFunc::ZEROS( porter.data(), GlobalC::pw.nrxx);
    for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
    {
        porter[GlobalC::pw.ig2fftc[ig]] = vg[ig];
    }
    GlobalC::pw.FFT_chg.FFT3D(porter.data(), 1);
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToReciSpace(const double* vr, std::complex<double> *vg)
{
	ModuleBase::GlobalFunc::ZEROS( porter.data(), GlobalC::pw.nrxx);
	for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		porter[ir] = std::complex<double>(vr[ir], 0.0);
	}
	GlobalC::pw.FFT_chg.FFT3D( porter.data(), -1);
	for (int ig=0; ig<GlobalC::pw.ngmc; ig++)
	{
		vg[ig] = porter[ GlobalC::pw.ig2fftc[ig] ];
	}
	return;
}
