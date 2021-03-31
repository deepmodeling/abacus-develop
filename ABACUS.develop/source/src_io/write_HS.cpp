#include "write_HS.h"
#include "../src_pw/global.h"


void HS_Matrix::saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_hs)
{   
    if(out_hs==1)
    {
        save_HS(Hloc, Sloc, bit);
    }
    else if(out_hs==2)
    {
        save_HS(Hloc, Sloc, bit);
    }
    else if(out_hs==3)
    {
        //please call individually
    }
    else if(out_hs==0)
    {
        // do nothing.
    }
    else
    {
        WARNING("Diago_LCAO_Matrix","unrecorganized out_hs value.");
    }
    return;
}


/*
void HS_Matrix::save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
        const double *nzval_H, const double *nzval_S, bool bit)
{
    TITLE("HS_Matrix","save_HS_ccf");

    if(DRANK!=0)return;
    
    stringstream ssh;
    stringstream sss;

    if(bit)
    {
        ssh << global_out_dir << "H_bit.ccf";
        sss << global_out_dir << "S_bit.ccf";
    }
    else
    {
		// mohan update 2021-02-10
        ssh << global_out_dir << "H" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
        sss << global_out_dir << "S" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
    }

    if(bit)
    {
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&NLOCAL,sizeof(int),1,g1);
        fwrite(&Hnnz,sizeof(int),1,g1);
        fwrite(&NLOCAL,sizeof(int),1,g2);
        fwrite(&Hnnz,sizeof(int),1,g2);

        fclose(g1);
        fclose(g2);
    }

        
    if(!bit)
    {
        ofstream g1(ssh.str().c_str());
        ofstream g2(sss.str().c_str());

        g1 << NLOCAL << " " << Hnnz << endl;
        g2 << NLOCAL << " " << Hnnz << endl;

        for(int i=0; i<NLOCAL+1; ++i)
        {
            g1 << colptr_H[i] << " ";
            g2 << colptr_H[i] << " ";
        }
        g1 << endl;
        g2 << endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << rowind_H[i] << " ";
            g2 << rowind_H[i] << " ";
        }
        g1 << endl;
        g2 << endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << nzval_H[i] << " ";
            g2 << nzval_S[i] << " ";
        }
        g1 << endl;
        g2 << endl;

        g1.close();
        g2.close();
    }
    return;
}
*/

// mohan add 2010/3/20, output H and S matrix, convinence for diagonalization
// test or save the middle information for next start.
void HS_Matrix::save_HS(const double *H, const double *S, bool bit)
{
    TITLE("HS_Matrix","save_HS_bit");
    timer::tick("HS_Matrix","save_HS_bit");
    OUT(ofs_running,"Dimension of H and S",NLOCAL);

    stringstream ssh;
    stringstream sss;

    if(bit)
    {
        ssh << global_out_dir << "data-H-bit";
        sss << global_out_dir << "data-S-bit";
    }
    else 
    {
        ssh << global_out_dir << "data-H";
        sss << global_out_dir << "data-S";
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (DRANK==0)
        {
            g1 = fopen(ssh.str().c_str(),"wb");
            g2 = fopen(sss.str().c_str(),"wb");
            fwrite(&NLOCAL,sizeof(int),1,g1);
            fwrite(&NLOCAL,sizeof(int),1,g2);
        }

        int ir,ic;
        for (int i=0; i<NLOCAL; i++)
        {
            double* lineH = new double[NLOCAL-i];
            double* lineS = new double[NLOCAL-i];
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
                        int iic;
                        if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*ParaO.nrow;
                        }
                        else
                        {
                            iic=ir*ParaO.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*ParaO.ncol+ic];
                        //lineS[j-i] = S[ir*ParaO.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,NLOCAL-i);

            if (DRANK==0)
            {
                for (int j=i; j<NLOCAL; j++)
                {
                    fwrite(&lineH[j-i],sizeof(double),1,g1);
                    fwrite(&lineS[j-i],sizeof(double),1,g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (DRANK==0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&NLOCAL,sizeof(int),1,g1);
        fwrite(&NLOCAL,sizeof(int),1,g2);

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                fwrite(&H[i*NLOCAL+j],sizeof(double),1,g1);
                fwrite(&S[i*NLOCAL+j],sizeof(double),1,g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } //end bit
    else
    {
#ifdef __MPI
        ofstream g1;
        ofstream g2;

        if (DRANK==0)
        {
            g1.open(ssh.str().c_str());
            g2.open(sss.str().c_str());
            g1 << NLOCAL;
            g2 << NLOCAL;
        }

        int ir,ic;
        for (int i=0; i<NLOCAL; i++)
        {
            double* lineH = new double[NLOCAL-i];
            double* lineS = new double[NLOCAL-i];
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
                        int iic;
                        if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*ParaO.nrow;
                        }
                        else
                        {
                            iic=ir*ParaO.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*ParaO.ncol+ic];
                        //lineS[j-i] = S[ir*ParaO.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,NLOCAL-i);

            if (DRANK==0)
            {
                for (int j=i; j<NLOCAL; j++)
                {
                    g1 << " " << lineH[j-i];
                    g2 << " " << lineS[j-i];
                }
                g1 << endl;
                g2 << endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        //if (DRANK==0);
        if (DRANK==0)       // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }

/*LiuXH add 2015-12-17,begin
    //int nprocs,myid;
    //MPI_Status status;
    //MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
    //MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

    string H_fn;
    stringstream H_fn2;
    H_fn2<< "data-H-"  << DRANK ;
    H_fn=H_fn2.str();
    ofstream ofs_H;
    ofs_H.open(H_fn.c_str());
    ofs_H<<setprecision(8) << setw(12);

    string S_fn;
    stringstream S_fn2;
    S_fn2<< "data-S-"  << DRANK ;
    S_fn=S_fn2.str();
    ofstream ofs_S;
    ofs_S.open(S_fn.c_str());
    ofs_S<<setprecision(8) << setw(12);

        int irr,icc;
        for (int i=0; i<NLOCAL; i++)
        {
            irr = ParaO.trace_loc_row[i];
            if (irr>=0)
            {
                // data collection
                for (int j=0; j<NLOCAL; j++)
                {
            icc = ParaO.trace_loc_col[j];
            if (icc>=0)
            {
                //if(abs(H[irr*ParaO.ncol+icc]) < 1.0e-10) H[irr*ParaO.ncol+icc] = 0.0;
                //if(abs(S[irr*ParaO.ncol+icc]) < 1.0e-10) S[irr*ParaO.ncol+icc] = 0.0;
                ofs_H << " " << H[irr*ParaO.ncol+icc];
                ofs_S << " " << S[irr*ParaO.ncol+icc];
            }
        }
        ofs_H << endl;
        ofs_S << endl;
         }
         }
//LiuXH add 2015-12-17,end*/
#else
        ofstream g1(ssh.str().c_str());
        ofstream g2(sss.str().c_str());

        g1 << NLOCAL;
        g2 << NLOCAL;

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                g1 << " " << H[i*NLOCAL+j];
                g2 << " " << S[i*NLOCAL+j];
            }
            g1 << endl;
            g2 << endl;
        }
        g1.close();
        g2.close();
#endif
    }

    timer::tick("HS_Matrix","save_HS_bit");
    return;
}

//LiuXh, 2017-03-21
void HS_Matrix::saving_HS_complex(complex<double> *Hloc, complex<double>* Sloc, bool bit, const int &out_hs)
{   
    if(out_hs==1)
    {
        save_HS_complex(Hloc, Sloc, bit);
    }
    else if(out_hs==0)
    {
        // do nothing.
    }
    else
    {
        WARNING("Diago_LCAO_Matrix","unrecorganized out_hs value.");
    }
    return;
}

//LiuXh, 2017-03-21
void HS_Matrix::save_HS_complex(complex<double> *H, complex<double> *S, bool bit)
{
    TITLE("HS_Matrix","save_HS_bit");
    timer::tick("HS_Matrix","save_HS_bit");
    OUT(ofs_running,"Dimension of H and S",NLOCAL);

    stringstream ssh;
    stringstream sss;

    if(bit)
    {
    ssh << global_out_dir << "data-H-bit";
    sss << global_out_dir << "data-S-bit";
    }
    else
    {
    ssh << global_out_dir << "data-H";
    sss << global_out_dir << "data-S";
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (DRANK==0)
        {
            g1 = fopen(ssh.str().c_str(),"wb");
            g2 = fopen(sss.str().c_str(),"wb");
            fwrite(&NLOCAL,sizeof(int),1,g1);
            fwrite(&NLOCAL,sizeof(int),1,g2);
        }

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
                        int iic;
                        if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*ParaO.nrow;
                        }
                        else
                        {
                            iic=ir*ParaO.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*ParaO.ncol+ic];
                        //lineS[j-i] = S[ir*ParaO.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
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
                    fwrite(&lineH[j-i],sizeof(complex<double>),1,g1);
                    fwrite(&lineS[j-i],sizeof(complex<double>),1,g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (DRANK==0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&NLOCAL,sizeof(int),1,g1);
        fwrite(&NLOCAL,sizeof(int),1,g2);

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                fwrite(&H[i*NLOCAL+j],sizeof(complex<double>),1,g1);
                fwrite(&S[i*NLOCAL+j],sizeof(complex<double>),1,g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } //end bit
    else
    {
#ifdef __MPI
        ofstream g1;
        ofstream g2;

        if (DRANK==0)
        {
            g1.open(ssh.str().c_str());
            g2.open(sss.str().c_str());
            g1 << NLOCAL;
            g2 << NLOCAL;
        }

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
                        int iic;
                        if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*ParaO.nrow;
                        }
                        else
                        {
                            iic=ir*ParaO.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*ParaO.ncol+ic];
                        //lineS[j-i] = S[ir*ParaO.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
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
                    g1 << " " << lineH[j-i];
                    g2 << " " << lineS[j-i];
                }
                g1 << endl;
                g2 << endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        //if (DRANK==0);
        if (DRANK==0)           // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }

/*LiuXH add 2015-12-17,begin
        //int nprocs,myid;
        //MPI_Status status;
        //MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
        //MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

        string H_fn;
        stringstream H_fn2;
        H_fn2<< "data-H-"  << DRANK ;
        H_fn=H_fn2.str();
        ofstream ofs_H;
        ofs_H.open(H_fn.c_str());
        ofs_H<<setprecision(8) << setw(12);

        string S_fn;
        stringstream S_fn2;
        S_fn2<< "data-S-"  << DRANK ;
        S_fn=S_fn2.str();
        ofstream ofs_S;
        ofs_S.open(S_fn.c_str());
        ofs_S<<setprecision(8) << setw(12);

        int irr,icc;
        for (int i=0; i<NLOCAL; i++)
        {
            irr = ParaO.trace_loc_row[i];
            if (irr>=0)
            {
                // data collection
                for (int j=0; j<NLOCAL; j++)
                {
                        icc = ParaO.trace_loc_col[j];
                        if (icc>=0)
                        {
                                //if(abs(H[irr*ParaO.ncol+icc]) < 1.0e-10) H[irr*ParaO.ncol+icc] = 0.0;
                                //if(abs(S[irr*ParaO.ncol+icc]) < 1.0e-10) S[irr*ParaO.ncol+icc] = 0.0;
                                ofs_H << " " << H[irr*ParaO.ncol+icc];
                                ofs_S << " " << S[irr*ParaO.ncol+icc];
                        }
                }
                ofs_H << endl;
                ofs_S << endl;
             }
         }
//LiuXH add 2015-12-17,end*/
#else
        ofstream g1(ssh.str().c_str());
        ofstream g2(sss.str().c_str());

        g1 << NLOCAL;
        g2 << NLOCAL;

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                g1 << " " << H[i*NLOCAL+j];
                g2 << " " << S[i*NLOCAL+j];
            }
            g1 << endl;
            g2 << endl;
        }
        g1.close();
        g2.close();
#endif
    }

    timer::tick("HS_Matrix","save_HS_bit");
    return;
}

//void HS_Matrix::save_HSR_tr(const int Rx, const int Ry, const int Rz, const double *H, const double *S)
void HS_Matrix::save_HSR_tr(const int current_spin)
//void HS_Matrix::save_HSR_tr(void)
{
    TITLE("HS_Matrix","save_HSR_tr");
    timer::tick("HS_Matrix","save_HSR_tr");

    stringstream ssh;
    stringstream sss;

    ssh << global_out_dir << "data-HR-tr_SPIN"<<current_spin;
    sss << global_out_dir << "data-SR-tr_SPIN"<<current_spin;
    //ssh << global_out_dir << "data-HR-tr_SPIN";
    //sss << global_out_dir << "data-SR-tr_SPIN";

#ifdef __MPI
    ofstream g1;
    ofstream g2;

    if(DRANK==0)
    {
        g1.open(ssh.str().c_str());
        g2.open(sss.str().c_str());
        g1 << "Matrix Dimension of H(R): "<<NLOCAL<<endl;
        g2 << "Matrix Dimension of S(R): "<<NLOCAL<<endl;
    }

    int R_x = GridD.getCellX();
    int R_y = GridD.getCellY();
    int R_z = GridD.getCellZ();

//cout<<"R_x: "<<R_x<<endl;
//cout<<"R_y: "<<R_y<<endl;
//cout<<"R_z: "<<R_z<<endl;

    double R_minX = GridD.getD_minX();
    double R_minY = GridD.getD_minY();
    double R_minZ = GridD.getD_minZ();

    //int dRx, dRy, dRz;

    for(int ix=0; ix<R_x; ix++)
    {
        int dRx = ix + R_minX;
        for(int iy=0; iy<R_y; iy++)
        {
            int dRy = iy + R_minY;
            for(int iz=0; iz<R_z; iz++)
            {
                int dRz = iz + R_minZ;
//cout<<"dRx: "<<dRx<<endl;
//cout<<"dRy: "<<dRy<<endl;
//cout<<"dRz: "<<dRz<<endl;
                int ir,ic;
                for(int i=0; i<NLOCAL; i++)
                {
                    //double* lineH = new double[NLOCAL-i];
                    //double* lineS = new double[NLOCAL-i];
                    double* lineH = nullptr;
                    double* lineS = nullptr;
                    complex<double>* lineH_soc = nullptr;
                    complex<double>* lineS_soc = nullptr;
                    if(NSPIN!=4)
                    {
                        lineH = new double[NLOCAL];
                        lineS = new double[NLOCAL];
                        ZEROS(lineH, NLOCAL);
                        ZEROS(lineS, NLOCAL);
                    }
                    else
                    {
                        lineH_soc = new complex<double>[NLOCAL];
                        lineS_soc = new complex<double>[NLOCAL];
                        ZEROS(lineH_soc, NLOCAL);
                        ZEROS(lineS_soc, NLOCAL);
                    }
                    //ZEROS(lineH, NLOCAL-i);
                    //ZEROS(lineS, NLOCAL-i);
                    //ZEROS(lineH, NLOCAL);
                    //ZEROS(lineS, NLOCAL);

                    ir = ParaO.trace_loc_row[i];
                    if(ir>=0)
                    {
                        //for(int j=i; j<NLOCAL; j++)
                        for(int j=0; j<NLOCAL; j++)
                        {
                            ic = ParaO.trace_loc_col[j];
                            if(ic>=0)
                            {
                                //lineH[j-i] = H[ir*ParaO.ncol+ic];
                                //lineS[j-i] = S[ir*ParaO.ncol+ic];
                                int iic;
                                if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                                {
                                    iic=ir+ic*ParaO.nrow;
                                }
                                else
                                {
                                    iic=ir*ParaO.ncol+ic;
                                }
                                if(NSPIN!=4)
                                {
                                    lineH[j] = LM.HR_tr[ix][iy][iz][iic];
                                    lineS[j] = LM.SlocR_tr[ix][iy][iz][iic];
                                }
                                else
                                {
                                    lineH_soc[j] = LM.HR_tr_soc[ix][iy][iz][iic];
                                    lineS_soc[j] = LM.SlocR_tr_soc[ix][iy][iz][iic];
                                }
                            }
                        }
                    }
                    else
                    {
                        //do nothing
                    }

                    //Parallel_Reduce::reduce_double_all(lineH,NLOCAL-i);
                    //Parallel_Reduce::reduce_double_all(lineS,NLOCAL-i);
                    if(NSPIN!=4)
                    {
                        Parallel_Reduce::reduce_double_all(lineH,NLOCAL);
                        Parallel_Reduce::reduce_double_all(lineS,NLOCAL);
                    }
                    else
                    {
                        Parallel_Reduce::reduce_complex_double_all(lineH_soc,NLOCAL);
                        Parallel_Reduce::reduce_complex_double_all(lineS_soc,NLOCAL);
                    }

                    if(DRANK==0)
                    {
                        //for(int j=i; j<NLOCAL; j++)
                        for(int j=0; j<NLOCAL; j++)
                        {
                            if(i==0 && j==0)
                            {
                                g1 << dRx << " " << dRy << " " << dRz  << "    //R vector(R2 - R1,unit: lattice vector)" <<endl;
                                g2 << dRx << " " << dRy << " " << dRz  << "    //R vector(R2 - R1,unit: lattice vector)" <<endl;
                            }
                            //g1 << " " << lineH[j-i];
                            //g2 << " " << lineS[j-i];
                            if(NSPIN!=4)
                            {
                                if(abs(lineH[j]) < 1.0e-12) lineH[j]=0.0;
                                if(abs(lineS[j]) < 1.0e-12) lineS[j]=0.0;
                                g1 << " " << lineH[j];
                                g2 << " " << lineS[j];
                            }
                            else
                            {
                if(abs(lineH_soc[j].real()) < 1.0e-12) lineH_soc[j]= complex<double> (0.0, lineH_soc[j].imag());
                                if(abs(lineH_soc[j].imag()) < 1.0e-12) lineH_soc[j]= complex<double> (lineH_soc[j].real(), 0.0);
                                if(abs(lineS_soc[j].real()) < 1.0e-12) lineS_soc[j]= complex<double> (0.0, lineS_soc[j].imag());
                                if(abs(lineS_soc[j].imag()) < 1.0e-12) lineS_soc[j]= complex<double> (lineS_soc[j].real() , 0.0);
                                g1 << " " << lineH_soc[j];
                                g2 << " " << lineS_soc[j];
                            }
                        }
                        g1 << endl;
                        g2 << endl;
                    }
                    if(NSPIN!=4)
                    {
                        delete[] lineH;
                        delete[] lineS;
                    }
                    else
                    {
                        delete[] lineH_soc;
                        delete[] lineS_soc;
                    }
                }
/*
                if(DRANK==0);
                {
                    g1.close();
                    g2.close();
                }
*/
            }
        }
    }
    //if(DRANK==0);
    if(DRANK==0)                // Peize Lin delete ; at 2020.01.31
    {
        g1.close();
        g2.close();
    }

#else
    ofstream g1(ssh.str().c_str());
    ofstream g2(sss.str().c_str());

    g1 << NLOCAL;
    g2 << NLOCAL;

    for (int i=0; i<NLOCAL; i++)
    {
        for (int j=i; j<NLOCAL; j++)
        {
            g1 << " " << H[i*NLOCAL+j];
            g2 << " " << S[i*NLOCAL+j];
        }
    }
    g1.close();
    g2.close();
#endif

    timer::tick("HS_Matrix","save_HSR_tr");
    return;
}
