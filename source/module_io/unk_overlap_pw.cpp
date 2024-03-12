#include "unk_overlap_pw.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"

unkOverlap_pw::unkOverlap_pw()
{
	//GlobalV::ofs_running << "this is unkOverlap_pw()" << std::endl;
}

unkOverlap_pw::~unkOverlap_pw()
{
	//GlobalV::ofs_running << "this is ~unkOverlap_pw()" << std::endl;
}

std::complex<double> unkOverlap_pw::unkdotp_G(const ModulePW::PW_Basis_K* wfcpw,
                                              const int ik_L,
                                              const int ik_R,
                                              const int iband_L,
                                              const int iband_R,
                                              const psi::Psi<std::complex<double>>* evc)
{
	
	std::complex<double> result(0.0,0.0);
	const int number_pw = wfcpw->npw;
	std::complex<double> *unk_L = new std::complex<double>[number_pw];
	std::complex<double> *unk_R = new std::complex<double>[number_pw];
	ModuleBase::GlobalFunc::ZEROS(unk_L,number_pw);
	ModuleBase::GlobalFunc::ZEROS(unk_R,number_pw);
	

	for (int igl = 0; igl < evc->get_ngk(ik_L); igl++)
	{
		unk_L[wfcpw->getigl2ig(ik_L,igl)] = evc[0](ik_L, iband_L, igl);
	}
	
	for (int igl = 0; igl < evc->get_ngk(ik_R); igl++)
	{
		unk_R[wfcpw->getigl2ig(ik_L,igl)] = evc[0](ik_R, iband_R, igl);
	}

	
	for (int iG = 0; iG < number_pw; iG++)
	{

		result = result + conj(unk_L[iG]) * unk_R[iG];

	}


#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = std::complex<double>(out_date_real,out_date_imag);
#endif

	delete[] unk_L;
	delete[] unk_R;
	return result;


}

std::complex<double> unkOverlap_pw::unkdotp_G0(const ModulePW::PW_Basis* rhopw,
                                               const ModulePW::PW_Basis_K* wfcpw,
                                               const int ik_L,
                                               const int ik_R,
                                               const int iband_L,
                                               const int iband_R,
                                               const psi::Psi<std::complex<double>>* evc,
                                               const ModuleBase::Vector3<double> G)
{
	// (1) set value
	std::complex<double> result(0.0,0.0);
    std::complex<double>* psi_r = new std::complex<double>[wfcpw->nmaxgr];
    std::complex<double>* phase = new std::complex<double>[rhopw->nmaxgr];

    // get the phase value in realspace
    for (int ig = 0; ig < rhopw->npw; ig++)
    {
		if (rhopw->gdirect[ig] == G)
		{
			phase[ig] = std::complex<double>(1.0,0.0);
			break;
		}
	}
	
	// (2) fft and get value
	rhopw->recip2real(phase, phase);
    wfcpw->recip2real(&evc[0](ik_L, iband_L, 0), psi_r, ik_L);

    for (int ir = 0; ir < rhopw->npw; ir++)
    {
        psi_r[ir] = psi_r[ir] * phase[ir];
    }

    // (3) calculate the overlap in ik_L and ik_R
    wfcpw->real2recip(psi_r, psi_r, ik_L);

    for (int ig = 0; ig < evc->get_ngk(ik_R); ig++)
    {
        result = result + conj(psi_r[ig]) * evc[0](ik_R, iband_R, ig);
    }

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = std::complex<double>(out_date_real,out_date_imag);
#endif
	
	delete[] psi_r;
	delete[] phase;
    return result;
}

// if noncollinear = 1 or GlobalV::NSPIN = 4 , you need this routine to calculate overlap unk
std::complex<double> unkOverlap_pw::unkdotp_soc_G(const ModulePW::PW_Basis_K* wfcpw,
                                                  const int ik_L,
                                                  const int ik_R,
                                                  const int iband_L,
                                                  const int iband_R,
                                                  const int npwx,
                                                  const psi::Psi<std::complex<double>>* evc)
{
	
	std::complex<double> result(0.0,0.0);
    const int number_pw = wfcpw->npw;
    std::complex<double>* unk_L = new std::complex<double>[number_pw * GlobalV::NPOL];
    std::complex<double>* unk_R = new std::complex<double>[number_pw * GlobalV::NPOL];
    ModuleBase::GlobalFunc::ZEROS(unk_L,number_pw*GlobalV::NPOL);
	ModuleBase::GlobalFunc::ZEROS(unk_R,number_pw*GlobalV::NPOL);
	
	for(int i = 0; i < GlobalV::NPOL; i++)
	{
		for (int igl = 0; igl < evc->get_ngk(ik_L); igl++)
		{
            unk_L[wfcpw->getigl2ig(ik_L, igl) + i * number_pw] = evc[0](ik_L, iband_L, igl + i * npwx);
        }

        for (int igl = 0; igl < evc->get_ngk(ik_R); igl++)
        {
            unk_R[wfcpw->getigl2ig(ik_L, igl) + i * number_pw] = evc[0](ik_R, iband_R, igl + i * npwx);
        }
    }

    for (int iG = 0; iG < number_pw*GlobalV::NPOL; iG++)
	{

		result = result + conj(unk_L[iG]) * unk_R[iG];

	}
	
#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = std::complex<double>(out_date_real,out_date_imag);
#endif

	delete[] unk_L;
	delete[] unk_R;
	return result;


}

//here G is in direct coordinate
std::complex<double> unkOverlap_pw::unkdotp_soc_G0(const ModulePW::PW_Basis* rhopw,
                                                   const ModulePW::PW_Basis_K* wfcpw,
                                                   const int ik_L,
                                                   const int ik_R,
                                                   const int iband_L,
                                                   const int iband_R,
                                                   const psi::Psi<std::complex<double>>* evc,
                                                   const ModuleBase::Vector3<double> G)
{
	// (1) set value
	std::complex<double> result(0.0,0.0);
	std::complex<double> *phase =new std::complex<double>[rhopw->nmaxgr];
    std::complex<double>* psi_up = new std::complex<double>[wfcpw->nmaxgr];
    std::complex<double>* psi_down = new std::complex<double>[wfcpw->nmaxgr];
    const int npwx = wfcpw->npwk_max;

    // get the phase value in realspace
    for (int ig = 0; ig < rhopw->npw; ig++)
	{
		if (rhopw->gdirect[ig] == G)
		{
			phase[ig] = std::complex<double>(1.0,0.0);
			break;
		}
	}
	
	// (2) fft and get value
	rhopw->recip2real(phase, phase);
    wfcpw->recip2real(&evc[0](ik_L, iband_L, 0), psi_up, ik_L);
    wfcpw->recip2real(&evc[0](ik_L, iband_L, npwx), psi_down, ik_L);

    for (int ir = 0; ir < wfcpw->nrxx; ir++)
    {
        psi_up[ir] = psi_up[ir] * phase[ir];
		psi_down[ir] = psi_down[ir] * phase[ir];
    }

    // (3) calculate the overlap in ik_L and ik_R
    wfcpw->real2recip(psi_up, psi_up, ik_L);
    wfcpw->real2recip(psi_down, psi_down, ik_L);

    for (int i = 0; i < GlobalV::NPOL; i++)
    {
		for(int ig = 0; ig < evc->get_ngk(ik_R); ig++)
		{
			if( i == 0 ) result = result + conj( psi_up[ig] ) * evc[0](ik_R, iband_R, ig);
			if( i == 1 ) result = result + conj( psi_down[ig] ) * evc[0](ik_R, iband_R, ig + npwx);
		}
	}
	
#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the GlobalV::KPAR = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = std::complex<double>(out_date_real,out_date_imag);
#endif
	
	delete[] psi_up;
	delete[] psi_down;
    return result;
}
