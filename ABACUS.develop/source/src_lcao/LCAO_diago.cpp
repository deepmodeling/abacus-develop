#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_io/write_HS.h"
//xiaohui modified 2013-03-23
//#include "../src_develop/src_siao/selinv.h"

Diago_LCAO_Matrix::Diago_LCAO_Matrix(){}
Diago_LCAO_Matrix::~Diago_LCAO_Matrix(){}

void Diago_LCAO_Matrix::solve_complex_matrix(
	const int &ik, 
	complex<double>** wfc, 
	ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_complex_matrix");
	time_t time_start = time(NULL);

	if(KS_SOLVER=="lapack")
	{
		this->using_LAPACK_complex(ik, wfc);
	}
	else
	{
#ifdef __MPI
		this->using_HPSEPS_complex(ik, wfc, wfc_2d);
#else
		WARNING_QUIT("Diago_LCAO_Matrix::solve_complex_matrix","only lapack is available!");
#endif
	}

	time_t time_end = time(NULL);

	OUT_TIME("diago(complex)", time_start, time_end);

	return;
}


void Diago_LCAO_Matrix::solve_double_matrix(
	const int &ik, 
	double** wfc, 
	matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","solve_double_matrix");
	timer::tick("Diago_LCAO_Matrix","solve_double_matrix",'F');
	time_t time_start = time(NULL);

	
	if(KS_SOLVER=="lapack")
	{
		this->using_LAPACK(ik, wfc);
	}
#ifdef __MPI
	else if(KS_SOLVER=="hpseps" || KS_SOLVER=="genelpa"|| KS_SOLVER=="scalapack_gvx")
	{
		this->using_HPSEPS_double(ik, wfc, wfc_2d);
	}
#endif
	else
	{
		cout << " Diago_LCAO_Matrix, diago_type = " << KS_SOLVER << endl; 
		WARNING_QUIT("Diago_LCAO_Matrix::init","Check KS_SOLVER.");
	}

	time_t time_end = time(NULL);

	OUT_TIME("diago(double)",time_start, time_end);

	timer::tick("Diago_LCAO_Matrix","solve_double_matrix",'F');
	return;
}

#ifdef __MPI
void Diago_LCAO_Matrix::using_HPSEPS_double(const int &ik, double**wfc, matrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_double");

	// save H and S matrix to disk.
	bool bit = false;
	HS_Matrix::saving_HS(LM.Hloc, LM.Sloc, bit, ParaO.out_hs);
	ofs_running << setprecision(6);

	// Distribution of matrix for 
	// prallel eigensolver.
	ParaO.diago_double_begin(ik, wfc, wfc_2d, LM.Hloc, LM.Sloc, wf.ekb[ik]);

	return;
}


void Diago_LCAO_Matrix::using_HPSEPS_complex(const int &ik, complex<double>** wfc, ComplexMatrix &wfc_2d)const
{
	TITLE("Diago_LCAO_Matrix","using_HPSEPS_complex");

	//ParaO.out_hs=1;//zhengdy-soc-test
	bool bit = false; //LiuXh, 2017-03-21
	//if set bit = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
	HS_Matrix::saving_HS_complex(LM.Hloc2, LM.Sloc2, bit, ParaO.out_hs); //LiuXh, 2017-03-21
	ofs_running << setprecision(6); //LiuXh, 2017-03-21

	ParaO.diago_complex_begin(ik, wfc, wfc_2d, LM.Hloc2, LM.Sloc2, wf.ekb[ik]);

	//added by zhengdy-soc, rearrange the WFC_K from [up,down,up,down...] to [up,up...down,down...], 
	if(NSPIN==4)
	{
		int row = GridT.lgd;
		vector<complex<double>> tmp(row);
		for(int ib=0; ib<NBANDS; ib++)
		{
			for(int iw=0; iw<row / NPOL; iw++)
			{
				tmp[iw] = LOWF.WFC_K[ik][ib][iw * NPOL];
				tmp[iw + row / NPOL] = LOWF.WFC_K[ik][ib][iw * NPOL + 1];
			}
			for(int iw=0; iw<row; iw++)
			{
				LOWF.WFC_K[ik][ib][iw] = tmp[iw];
			}
		}
	}

	return;
}
#endif

void Diago_LCAO_Matrix::using_LAPACK_complex(const int &ik, complex<double> **wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK_complex");

	assert(NPROC = 1);

	ComplexMatrix Htmp(NLOCAL,NLOCAL);
	ComplexMatrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			Htmp(i,j) = LM.Hloc2[i*NLOCAL+j];
			Stmp(i,j) = LM.Sloc2[i*NLOCAL+j];
		}
	}

//----------------------------
// keep this for tests
//    out.printcm_norm("Lapack_H", Htmp, 1.0e-5);
//    out.printcm_norm("Lapack_S", Stmp, 1.0e-5);
//----------------------------

	double* en = new double[NLOCAL];
	ZEROS(en, NLOCAL);

	ComplexMatrix hvec(NLOCAL, NBANDS);
	hm.cdiaghg(NLOCAL, NBANDS, Htmp, Stmp, NLOCAL, en, hvec);

	if(NSPIN!=4)
	{
		for(int ib=0; ib<NBANDS; ib++)
		{
			for(int iw=0; iw<NLOCAL; iw++)
			{
				LOWF.WFC_K[ik][ib][iw] = hvec(iw,ib);
			}
		}
	}
	else
	{
		for(int ib=0; ib<NBANDS; ib++)
		{
			for(int iw=0; iw<NLOCAL / NPOL; iw++)
			{
				LOWF.WFC_K[ik][ib][iw] = hvec(iw * NPOL, ib);
				LOWF.WFC_K[ik][ib][iw + NLOCAL / NPOL] = hvec(iw * NPOL + 1, ib);
			}
		}
	}

	// energy for k-point ik
	for(int ib=0; ib<NBANDS; ib++)
	{
		wf.ekb[ik][ib] = en[ib];
	}

	return;
}


void Diago_LCAO_Matrix::using_LAPACK(const int &ik, double** wfc)const
{
	TITLE("Diago_LCAO_Matrix","using_LAPACK");
	assert(NLOCAL>0);

	// save H and S matrix to disk.
//	bool bit = false;
	bool bit = true;//zhengdy-soc
	HS_Matrix::saving_HS(LM.Hloc, LM.Sloc, bit, ParaO.out_hs);

	matrix Htmp(NLOCAL,NLOCAL);
	matrix Stmp(NLOCAL,NLOCAL);
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			Htmp(i,j) = LM.Hloc[i*NLOCAL+j];
			Stmp(i,j) = LM.Sloc[i*NLOCAL+j];
		}
	}

	// @@@@@@@
	// test
	// @@@@@@@
//    out.printrm("Lapack_H", Htmp);
//    out.printrm("Lapack_S", Stmp);
		
	int itype=1;
	int lwork=3*NLOCAL-1;// tmp
	double* w = new double[NLOCAL];
	double* work = new double[lwork];
	ZEROS(w, NLOCAL);
	ZEROS(work, lwork);
	int info;

	clock_t clock_start, clock_end;
	clock_start = std::clock();
	LapackConnector::dsygv(itype,'V','U',NLOCAL,Htmp,NLOCAL,Stmp,NLOCAL,w,work,lwork,&info);
	clock_end = std::clock();
	double duration = (double)(clock_end - clock_start) / CLOCKS_PER_SEC;

	ofs_running << setiosflags(ios::fixed) << setprecision(20);
//	ofs_running << " clock_start = " << clock_start << endl;
//	ofs_running << " clock_end = " << clock_end << endl;
	ofs_running << " Time using dsygv in LAPACK (seconds) is " << duration << endl;
	ofs_running << resetiosflags(ios::fixed) << setprecision(10);

	for(int i=0; i<NBANDS; i++)
	{
		// eigenvalues
		wf.ekb[ik][i] = w[i]; 
		for(int j=0; j<NLOCAL; j++)
		{
			wfc[i][j] = Htmp(j,i);
		}
	}
	

	// @@@@@@@
	// test
	// @@@@@@@
	/*
	cout << "\n Lapack, wfc after diago:" << endl;
	for(int i=0; i<NBANDS; i++)
	{
		cout << " Eigenvalue from LAPACK : " << setw(5) << setw(12) << wf.ekb[ik][i] << endl;
		cout << " Eigenfunctions" << endl;
		for(int j=0; j<NLOCAL; j++)
		{
			cout << setw(12) << wfc[i][j];
		}
		cout << endl;
	}
	//exit(0);
	*/

	delete[] w;
	delete[] work;
	return;
}
