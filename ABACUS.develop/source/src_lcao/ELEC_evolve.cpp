#include "ELEC_evolve.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "src_pw/global.h"
#include "src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"

ELEC_evolve::ELEC_evolve(){};
ELEC_evolve::~ELEC_evolve(){};

// this routine only serves for TDDFT using LCAO basis set
void ELEC_evolve::evolve_psi(const int &istep, LCAO_Hamilt &uhm, complex<double> ***wfc)
{
	TITLE("ELEC_evolve","eveolve_psi");
	timer::tick("ELEC_evolve","evolve_psi",'E');

	int start_spin = -1;
	uhm.GK.reset_spin(start_spin);
	uhm.GK.allocate_pvpR();
						
	// pool parallization in future -- mohan note 2021-02-09
	for(int ik=0; ik<kv.nks; ik++)
	{	
		//-----------------------------------------
		//(1) prepare data for this k point.
		// copy the local potential from array.
		//-----------------------------------------
		if(NSPIN==2) 
		{
			CURRENT_SPIN = kv.isk[ik];
		}
		wf.npw = kv.ngk[ik];
		for(int ir=0; ir<pw.nrxx; ir++)
		{
			pot.vr_eff1[ir] = pot.vr_eff( CURRENT_SPIN, ir);
		}
		
		//--------------------------------------------
		//(2) check if we need to calculate 
		// pvpR = < phi0 | v(spin) | phiR> for a new spin.
		//--------------------------------------------
		if(CURRENT_SPIN == uhm.GK.get_spin() )
		{
			//ofs_running << " Same spin, same vlocal integration." << endl;
		}
		else
		{
			uhm.GK.reset_spin( CURRENT_SPIN );

			// vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
			uhm.GK.cal_vlocal_k(pot.vr_eff1,GridT);
			// added by zhengdy-soc, for non-collinear case
			// integral 4 times, is there any method to simplify?
			if(NSPIN==4)
			{
				for(int is=1;is<4;is++)
				{
					for(int ir=0; ir<pw.nrxx; ir++)
					{
						pot.vr_eff1[ir] = pot.vr_eff(is, ir);
					}
					uhm.GK.cal_vlocal_k(pot.vr_eff1, GridT, is);
				}
			}
		}


		if(!uhm.init_s)
    	{
    	    WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// with k points
		uhm.calculate_Hk(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u)
		{		
			dftu.cal_eff_pot_mat(ik, istep);

			for(int irc=0; irc<ParaO.nloc; irc++)
			{
				LM.Hloc2[irc] += dftu.pot_eff_k.at(ik).at(irc);
			}							
		}

		// Peize Lin add at 2020.04.04
		if(restart.info_load.load_H && !restart.info_load.load_H_finish)
		{
			restart.load_disk("H", ik);
			restart.info_load.load_H_finish = true;
		}			
		if(restart.info_save.save_H)
		{
			restart.save_disk("H", ik);
		}

		bool diago = true;
		if (istep >= 1) diago = false;
		if(diago)
		{
			timer::tick("Efficience","diago_k");
			Diago_LCAO_Matrix DLM;
			DLM.solve_complex_matrix(ik, LOWF.WFC_K[ik], LOC.wfc_dm_2d.wfc_k[ik]);
			timer::tick("Efficience","diago_k");
		}
		else
		{
			timer::tick("Efficience","evolve_k");
			Evolve_LCAO_Matrix ELM;
			ELM.evolve_complex_matrix(ik, LOWF.WFC_K[ik], wfc[ik]);
			timer::tick("Efficience","evolve_k");
		}
	} // end k
			
	// LiuXh modify 2019-07-15*/
	if(!ParaO.out_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	timer::tick("ELEC_evolve","evolve_psi",'E');
	return;	
}


void ELEC_evolve::evolve_complex_matrix(const int &ik, complex<double>** cc, complex<double>** cc_init)const
{
	TITLE("Evolve_LCAO_Matrix","evolve_complex_matrix");
	time_t time_start = time(NULL);
	ofs_running << " Start Time : " << ctime(&time_start);

	if (tddft==1)
	{
/*
#ifdef __MPI
		this->using_LAPACK_complex_2(ik, cc, cc_init);
#else
		this->using_LAPACK_complex(ik, cc, cc_init);
#endif
*/
		this-> using_LAPACK_complex(ik, cc, cc_init);
	}
	else
	{
		WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix","only tddft==1 cando evolve");
	}

	time_t time_end = time(NULL);
	OUT_TIME("evolve(complex)", time_start, time_end);
	
	return;
}

void ELEC_evolve::using_LAPACK_complex(const int &ik, complex<double>** c, complex<double>** c_init)const
{
	                                                                                                                     TITLE("Evolve_LCAO_Matrix","using_LAPACK_complex");

//	Calculate the U operator

	bool bit = false;
	//HS_Matrix::saving_HS_complex(LM.Hloc2, LM.Sloc2, bit, ParaO.out_hs);

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

				
	int INFO;

        int LWORK=3*NLOCAL-1; //tmp
        complex<double> * WORK = new complex<double>[LWORK];
        ZEROS(WORK, LWORK);
        int IPIV[NLOCAL];

        LapackConnector::zgetrf( NLOCAL, NLOCAL, Stmp, NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri( NLOCAL, Stmp, NLOCAL, IPIV, WORK, LWORK, &INFO);


	ComplexMatrix S_plus_H(NLOCAL,NLOCAL);
	S_plus_H = Stmp*Htmp;


	ComplexMatrix Denominator(NLOCAL,NLOCAL);
	for (int i=0; i<NLOCAL; i++)
       	{
               	for (int j=0; j<NLOCAL; j++)
                {
                     /*   real(Denominator(i,j)) = -imag(S_plus_H(i,j));
                        imag(Denominator(i,j)) = real(S_plus_H(i,j));*/
                          Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                }
        }
        
        ComplexMatrix Idmat(NLOCAL,NLOCAL);
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        if(i==j) Idmat(i,j) = complex<double>(1.0, 0.0);
                       	else Idmat(i,j) = complex<double>(0.0, 0.0);
                }
        }
        double delta_t;
//      delta_t = 0.2;	//identity: fs;
        ComplexMatrix Numerator(NLOCAL,NLOCAL);
        Numerator = Idmat - 0.5*INPUT.md_dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.md_dt*41.34*Denominator;

	int info;
        int lwork=3*NLOCAL-1; //tmp
        complex<double> * work = new complex<double>[lwork];
        ZEROS(work, lwork);
        int ipiv[NLOCAL];
        
        LapackConnector::zgetrf( NLOCAL, NLOCAL, Denominator, NLOCAL, ipiv, &info);
        LapackConnector::zgetri( NLOCAL, Denominator, NLOCAL, ipiv, work, lwork, &info);

        ComplexMatrix U_operator(NLOCAL,NLOCAL);

        U_operator = Numerator*Denominator;

// Calculate wave function at t+delta t
				
//	cout << "wave function coe at t+delta t !" << endl;

/*	for(int i=0; i<NBANDS; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			cout << c[i][j] << "\t";
		}
		cout <<endl;
	}
	cout << endl;
*/

	for(int i=0; i<NBANDS; i++)
	{
		complex<double> ccc[NLOCAL];
		for(int j=0; j<NLOCAL; j++)
		{	
			ccc[j] = (0.0,0.0);
			for(int k=0; k<NLOCAL; k++)
			{
				 ccc[j] += U_operator(j,k)*c_init[i][k];
			}
		}
		for(int j=0; j<NLOCAL; j++)
		{
			c[i][j] = ccc[j];
			LOWF.WFC_K[ik][i][j] = ccc[j];
		}	
	}

/*	for(int i=0; i<NBANDS; i++)
	{
                for(int j=0; j<NLOCAL; j++)
		{
			cout << LOWF.WFC_K[ik][i][j] << "\t";
		}
		cout <<endl;
	}
*/

//	delete[] work;
//	delete[] ipiv;

	return;
}

void ELEC_evolve::using_LAPACK_complex_2(const int &ik, complex<double>** c, complex<double>** c_init)const
{
	ComplexMatrix Htmp(NLOCAL,NLOCAL);
	ComplexMatrix Stmp(NLOCAL,NLOCAL);

        int ir,ic;
        for (int i=0; i<NLOCAL; i++)
        {
            complex<double>* lineH = new complex<double>[NLOCAL-i];
            complex<double>* lineS = new complex<double>[NLOCAL-i];
            ZEROS(lineH, NLOCAL-i);
            ZEROS(lineS, NLOCAL-i);

            ir = ParaO.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<NLOCAL; j++)
                {
                    ic = ParaO.trace_loc_col[j];
                    if (ic>=0)
                    {
                        //lineH[j-i] = H[ir*ParaO.ncol+ic];
			lineH[j-i] = LM.Hloc2[ir*ParaO.ncol+ic];
                        //lineS[j-i] = S[ir*ParaO.ncol+ic];
			lineS[j-i] = LM.Sloc2[ir*ParaO.ncol+ic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH,NLOCAL-i);
            Parallel_Reduce::reduce_complex_double_pool(lineS,NLOCAL-i);

            if (DRANK==0)
            {
                for (int j=i; j<NLOCAL; j++)
                {
                    //g1 << " " << lineH[j-i];
			Htmp(i,j) = lineH[j-i];
                    //g2 << " " << lineS[j-i];
			Stmp(i,j) = lineS[j-i];
                }
            }
            delete[] lineH;
            delete[] lineS;
        }
	
	int INFO;

        int LWORK=3*NLOCAL-1; //tmp
        complex<double> * WORK = new complex<double>[LWORK];
        ZEROS(WORK, LWORK);
        int IPIV[NLOCAL];

        LapackConnector::zgetrf( NLOCAL, NLOCAL, Stmp, NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri( NLOCAL, Stmp, NLOCAL, IPIV, WORK, LWORK, &INFO);

	ComplexMatrix S_plus_H(NLOCAL,NLOCAL);
	S_plus_H = Stmp*Htmp;

	ComplexMatrix Denominator(NLOCAL,NLOCAL);
	for (int i=0; i<NLOCAL; i++)
       	{
               	for (int j=0; j<NLOCAL; j++)
                {
                     /*   real(Denominator(i,j)) = -imag(S_plus_H(i,j));
                        imag(Denominator(i,j)) = real(S_plus_H(i,j));*/
                          Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                }
        }
        
        ComplexMatrix Idmat(NLOCAL,NLOCAL);
        for(int i=0; i<NLOCAL; i++)
        {
                for(int j=0; j<NLOCAL; j++)
                {
                        if(i==j) Idmat(i,j) = complex<double>(1.0, 0.0);
                       	else Idmat(i,j) = complex<double>(0.0, 0.0);
                }
        }
        double delta_t;
        ComplexMatrix Numerator(NLOCAL,NLOCAL);
        Numerator = Idmat - 0.5*INPUT.md_dt*41.34*Denominator;
        Denominator = Idmat + 0.5*INPUT.md_dt*41.34*Denominator;

	int info;
        int lwork=3*NLOCAL-1; //tmp
        complex<double> * work = new complex<double>[lwork];
        ZEROS(work, lwork);
        int ipiv[NLOCAL];
        
        LapackConnector::zgetrf( NLOCAL, NLOCAL, Denominator, NLOCAL, ipiv, &info);
        LapackConnector::zgetri( NLOCAL, Denominator, NLOCAL, ipiv, work, lwork, &info);

        ComplexMatrix U_operator(NLOCAL,NLOCAL);
        U_operator = Numerator*Denominator;

// Calculate wave function at t+delta t
				
//	cout << "wave function coe at t+delta t !" << endl;

	for(int i=0; i<NBANDS; i++)
	{
		complex<double> ccc[NLOCAL];
		for(int j=0; j<NLOCAL; j++)
		{	
			ccc[j] = (0.0,0.0);
			for(int k=0; k<NLOCAL; k++)
			{
				 ccc[j] += U_operator(j,k)*c_init[i][k];
			}
		}
		for(int j=0; j<NLOCAL; j++)
		{
			c[i][j] = ccc[j];
			LOWF.WFC_K[ik][i][j] = ccc[j];
		}	
	}

	return;
}
