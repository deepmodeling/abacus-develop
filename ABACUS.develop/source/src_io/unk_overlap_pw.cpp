#include "unk_overlap_pw.h"

unkOverlap_pw::unkOverlap_pw()
{
	//ofs_running << "this is unkOverlap_pw()" << endl;
}

unkOverlap_pw::~unkOverlap_pw()
{
	//ofs_running << "this is ~unkOverlap_pw()" << endl;
}



complex<double> unkOverlap_pw::unkdotp_G(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ComplexMatrix *evc)
{
	
	complex<double> result(0.0,0.0);
	// Total number of plane wave basis sets of wave function
	const int number_pw = pw.ngmw;
	complex<double> *unk_L = new complex<double>[number_pw];
	complex<double> *unk_R = new complex<double>[number_pw];
	ZEROS(unk_L,number_pw);
	ZEROS(unk_R,number_pw);
	

	for (int ig = 0; ig < kv.ngk[ik_L]; ig++)
	{
		unk_L[wf.igk(ik_L,ig)] = evc[ik_L](iband_L, ig);
	}
	
	for (int ig = 0; ig < kv.ngk[ik_R]; ig++)
	{
		unk_R[wf.igk(ik_R,ig)] = evc[ik_R](iband_R, ig);
	}

	
	for (int iG = 0; iG < number_pw; iG++)
	{

		result = result + conj(unk_L[iG]) * unk_R[iG];

	}


#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = complex<double>(out_date_real,out_date_imag);
#endif

	delete[] unk_L;
	delete[] unk_R;
	return result;


}



complex<double> unkOverlap_pw::unkdotp_G0(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ComplexMatrix *evc, const Vector3<double> G)
{
	// (1) set value
	complex<double> result(0.0,0.0);
	complex<double> *phase = UFFT.porter;
	complex<double> *psi_r = new complex<double>[pw.nrxx]; // Wave function in real space

	ZEROS( phase, pw.nrxx);
	ZEROS( psi_r, pw.nrxx );


	for (int ig = 0; ig < kv.ngk[ik_L]; ig++)
	{
		psi_r[ pw.ig2fftw[ wf.igk(ik_L,ig) ] ] = evc[ik_L](iband_L, ig);

	}

	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psi_r, 1);
	pw.FFT_wfc.FFT3D(phase, 1);
	
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psi_r[ir] = psi_r[ir] * phase[ir];
	}
	
	
	// (3) calculate the overlap in ik_L and ik_R
	pw.FFT_wfc.FFT3D( psi_r, -1);
	

	for(int ig = 0; ig < kv.ngk[ik_R]; ig++)
	{
		result = result + conj( psi_r[ pw.ig2fftw[wf.igk(ik_R,ig)] ] ) * evc[ik_R](iband_R, ig);
	}

#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = complex<double>(out_date_real,out_date_imag);
#endif
	
	delete[] psi_r;
    return result;
}

// if noncollinear = 1 or NSPIN = 4 , you need this routine to calculate overlap unk
complex<double> unkOverlap_pw::unkdotp_soc_G(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ComplexMatrix *evc)
{
	
	complex<double> result(0.0,0.0);
	//Total number of plane wave basis sets of wave function
	const int number_pw = pw.ngmw;
	complex<double> *unk_L = new complex<double>[number_pw*NPOL];
	complex<double> *unk_R = new complex<double>[number_pw*NPOL];
	ZEROS(unk_L,number_pw*NPOL);
	ZEROS(unk_R,number_pw*NPOL);
	
	for(int i = 0; i < NPOL; i++)
	{
		for (int ig = 0; ig < kv.ngk[ik_L]; ig++)
		{
			unk_L[wf.igk(ik_L,ig)+i*number_pw] = evc[ik_L](iband_L, ig+i*wf.npwx);
		}
	
		for (int ig = 0; ig < kv.ngk[ik_R]; ig++)
		{
			unk_R[wf.igk(ik_R,ig)+i*number_pw] = evc[ik_R](iband_R, ig+i*wf.npwx);
		}

	}
	
	for (int iG = 0; iG < number_pw*NPOL; iG++)
	{

		result = result + conj(unk_L[iG]) * unk_R[iG];

	}
	
#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = complex<double>(out_date_real,out_date_imag);
#endif

	delete[] unk_L;
	delete[] unk_R;
	return result;


}


// Here the G vector is the direct coordinate
complex<double> unkOverlap_pw::unkdotp_soc_G0(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ComplexMatrix *evc, const Vector3<double> G)
{
	// (1) set value
	complex<double> result(0.0,0.0);
	complex<double> *phase = UFFT.porter;
	complex<double> *psi_up = new complex<double>[pw.nrxx];
	complex<double> *psi_down = new complex<double>[pw.nrxx];
	ZEROS( phase, pw.nrxx);
	ZEROS( psi_up, pw.nrxx );
	ZEROS( psi_down, pw.nrxx );

	for(int i = 0; i < NPOL; i++)
	{
		for (int ig = 0; ig < kv.ngk[ik_L]; ig++)
		{
			if( i == 0 ) psi_up[ pw.ig2fftw[ wf.igk(ik_L,ig) ] ] = evc[ik_L](iband_L, ig);
			if( i == 1 ) psi_down[ pw.ig2fftw[ wf.igk(ik_L,ig) ] ] = evc[ik_L](iband_L, ig+i*wf.npwx);
		}
	}
	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psi_up, 1);
    pw.FFT_wfc.FFT3D(psi_down, 1);	
	pw.FFT_wfc.FFT3D(phase, 1);
	
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psi_up[ir] = psi_up[ir] * phase[ir];
		psi_down[ir] = psi_down[ir] * phase[ir];
	}
	
	
	// (3) calculate the overlap in ik_L and ik_R
	pw.FFT_wfc.FFT3D( psi_up, -1);
	pw.FFT_wfc.FFT3D( psi_down, -1);
	
	for(int i = 0; i < NPOL; i++)
	{
		for(int ig = 0; ig < kv.ngk[ik_R]; ig++)
		{
			if( i == 0 ) result = result + conj( psi_up[ pw.ig2fftw[wf.igk(ik_R,ig)] ] ) * evc[ik_R](iband_R, ig);
			if( i == 1 ) result = result + conj( psi_down[ pw.ig2fftw[wf.igk(ik_R,ig)] ] ) * evc[ik_R](iband_R, ig+i*wf.npwx);
		}
	}
	
#ifdef __MPI
    // note: the mpi uses MPI_COMMON_WORLD,so you must make the NPOOL = 1.
	double in_date_real = result.real();
	double in_date_imag = result.imag();
	double out_date_real = 0.0;
	double out_date_imag = 0.0;
	MPI_Allreduce(&in_date_real , &out_date_real , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	MPI_Allreduce(&in_date_imag , &out_date_imag , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
	result = complex<double>(out_date_real,out_date_imag);
#endif
	
	delete[] psi_up;
	delete[] psi_down;
    return result;
}

void unkOverlap_pw::test_for_unkOverlap_pw()
{
	
	const int number_pw = pw.ngmw;
	ofs_running << "the pw.ngmw is " << number_pw << endl;
	complex<double> *unk_L = new complex<double>[number_pw];
	for (int ig = 0; ig < kv.ngk[0]; ig++)
	{
		unk_L[wf.igk(0,ig)] = wf.evc[0](0, ig);
	}
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		ofs_running << pw.gdirect[ig].x << "," << pw.gdirect[ig].y << "," << pw.gdirect[ig].z << "  = " << unk_L[ig] << endl;
	}	
	
}

