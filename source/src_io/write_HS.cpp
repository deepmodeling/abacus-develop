#include "write_HS.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

void HS_Matrix::saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_mat_hs,  const Parallel_Orbitals &pv)
{   
    if(out_mat_hs==1)
    {
        save_HS(Hloc, Sloc, bit, pv);
    }
    else if(out_mat_hs==2)
    {
        save_HS(Hloc, Sloc, bit, pv);
    }
    else if(out_mat_hs==3)
    {
        //please call individually
    }
    else if(out_mat_hs==0)
    {
        // do nothing.
    }
    else
    {
        ModuleBase::WARNING("Diago_LCAO_Matrix","unrecorganized out_mat_hs value.");
    }
    return;
}


/*
void HS_Matrix::save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
        const double *nzval_H, const double *nzval_S, bool bit)
{
    ModuleBase::TITLE("HS_Matrix","save_HS_ccf");

    if(GlobalV::DRANK!=0)return;
    
    std::stringstream ssh;
    std::stringstream sss;

    if(bit)
    {
        ssh << GlobalV::global_out_dir << "H_bit.ccf";
        sss << GlobalV::global_out_dir << "S_bit.ccf";
    }
    else
    {
        // mohan update 2021-02-10
        ssh << GlobalV::global_out_dir << "H" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
        sss << GlobalV::global_out_dir << "S" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
    }

    if(bit)
    {
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
        fwrite(&Hnnz,sizeof(int),1,g1);
        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);
        fwrite(&Hnnz,sizeof(int),1,g2);

        fclose(g1);
        fclose(g2);
    }

        
    if(!bit)
    {
        std::ofstream g1(ssh.str().c_str());
        std::ofstream g2(sss.str().c_str());

        g1 << GlobalV::NLOCAL << " " << Hnnz << std::endl;
        g2 << GlobalV::NLOCAL << " " << Hnnz << std::endl;

        for(int i=0; i<GlobalV::NLOCAL+1; ++i)
        {
            g1 << colptr_H[i] << " ";
            g2 << colptr_H[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << rowind_H[i] << " ";
            g2 << rowind_H[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << nzval_H[i] << " ";
            g2 << nzval_S[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        g1.close();
        g2.close();
    }
    return;
}
*/

// mohan add 2010/3/20, output H and S matrix, convinence for diagonalization
// test or save the middle information for next start.
void HS_Matrix::save_HS(const double *H, const double *S, bool bit, const Parallel_Orbitals &pv)
{
    ModuleBase::TITLE("HS_Matrix","save_HS_bit");
    ModuleBase::timer::tick("HS_Matrix","save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Dimension of H and S",GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if(bit)
    {
        ssh << GlobalV::global_out_dir << "data-H-bit";
        sss << GlobalV::global_out_dir << "data-S-bit";
    }
    else 
    {
        ssh << GlobalV::global_out_dir << "data-H";
        sss << GlobalV::global_out_dir << "data-S";
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK==0)
        {
            g1 = fopen(ssh.str().c_str(),"wb");
            g2 = fopen(sss.str().c_str(),"wb");
            fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
            fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);
        }

        int ir,ic;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            double* lineH = new double[GlobalV::NLOCAL-i];
            double* lineS = new double[GlobalV::NLOCAL-i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);

            ir = pv.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    ic = pv.trace_loc_col[j];
                    if (ic>=0)
                    {
                        int iic;
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*pv.nrow;
                        }
                        else
                        {
                            iic=ir*pv.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*pv.ncol+ic];
                        //lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,GlobalV::NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,GlobalV::NLOCAL-i);

            if (GlobalV::DRANK==0)
            {
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j-i],sizeof(double),1,g1);
                    fwrite(&lineS[j-i],sizeof(double),1,g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK==0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            for (int j=i; j<GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i*GlobalV::NLOCAL+j],sizeof(double),1,g1);
                fwrite(&S[i*GlobalV::NLOCAL+j],sizeof(double),1,g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } //end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK==0)
        {
            g1.open(ssh.str().c_str());
            g2.open(sss.str().c_str());
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir,ic;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            double* lineH = new double[GlobalV::NLOCAL-i];
            double* lineS = new double[GlobalV::NLOCAL-i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);

            ir = pv.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    ic = pv.trace_loc_col[j];
                    if (ic>=0)
                    {
                        int iic;
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*pv.nrow;
                        }
                        else
                        {
                            iic=ir*pv.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*pv.ncol+ic];
                        //lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,GlobalV::NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,GlobalV::NLOCAL-i);

            if (GlobalV::DRANK==0)
            {
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j-i];
                    g2 << " " << lineS[j-i];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        //if (GlobalV::DRANK==0);
        if (GlobalV::DRANK==0)       // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }

/*LiuXH add 2015-12-17,begin
    //int nprocs,myid;
    //MPI_Status status;
    //MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
    //MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

    std::string H_fn;
    std::stringstream H_fn2;
    H_fn2<< "data-H-"  << GlobalV::DRANK ;
    H_fn=H_fn2.str();
    std::ofstream ofs_H;
    ofs_H.open(H_fn.c_str());
    ofs_H<<std::setprecision(8) << std::setw(12);

    std::string S_fn;
    std::stringstream S_fn2;
    S_fn2<< "data-S-"  << GlobalV::DRANK ;
    S_fn=S_fn2.str();
    std::ofstream ofs_S;
    ofs_S.open(S_fn.c_str());
    ofs_S<<std::setprecision(8) << std::setw(12);

        int irr,icc;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            irr = pv.trace_loc_row[i];
            if (irr>=0)
            {
                // data collection
                for (int j=0; j<GlobalV::NLOCAL; j++)
                {
            icc = pv.trace_loc_col[j];
            if (icc>=0)
            {
                //if(abs(H[irr*pv.ncol+icc]) < 1.0e-10) H[irr*pv.ncol+icc] = 0.0;
                //if(abs(S[irr*pv.ncol+icc]) < 1.0e-10) S[irr*pv.ncol+icc] = 0.0;
                ofs_H << " " << H[irr*pv.ncol+icc];
                ofs_S << " " << S[irr*pv.ncol+icc];
            }
        }
        ofs_H << std::endl;
        ofs_S << std::endl;
         }
         }
//LiuXH add 2015-12-17,end*/
#else
        std::ofstream g1(ssh.str().c_str());
        std::ofstream g2(sss.str().c_str());

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            for (int j=i; j<GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i*GlobalV::NLOCAL+j];
                g2 << " " << S[i*GlobalV::NLOCAL+j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("HS_Matrix","save_HS_bit");
    return;
}

//LiuXh, 2017-03-21
void HS_Matrix::saving_HS_complex(std::complex<double> *Hloc, std::complex<double>* Sloc, bool bit, const int &out_mat_hs, const Parallel_Orbitals &pv)
{   
    if(out_mat_hs==1)
    {
        save_HS_complex(Hloc, Sloc, bit, pv);
    }
    else if(out_mat_hs==0)
    {
        // do nothing.
    }
    else
    {
        ModuleBase::WARNING("Diago_LCAO_Matrix","unrecorganized out_mat_hs value.");
    }
    return;
}

//LiuXh, 2017-03-21
void HS_Matrix::save_HS_complex(std::complex<double> *H, std::complex<double> *S, bool bit, const Parallel_Orbitals &pv)
{
    ModuleBase::TITLE("HS_Matrix","save_HS_bit");
    ModuleBase::timer::tick("HS_Matrix","save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Dimension of H and S",GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if(bit)
    {
    ssh << GlobalV::global_out_dir << "data-H-bit";
    sss << GlobalV::global_out_dir << "data-S-bit";
    }
    else
    {
    ssh << GlobalV::global_out_dir << "data-H";
    sss << GlobalV::global_out_dir << "data-S";
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK==0)
        {
            g1 = fopen(ssh.str().c_str(),"wb");
            g2 = fopen(sss.str().c_str(),"wb");
            fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
            fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);
        }

        int ir,ic;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            std::complex<double>* lineH = new std::complex<double>[GlobalV::NLOCAL-i];
            std::complex<double>* lineS = new std::complex<double>[GlobalV::NLOCAL-i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);

            ir = pv.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    ic = pv.trace_loc_col[j];
                    if (ic>=0)
                    {
                        int iic;
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*pv.nrow;
                        }
                        else
                        {
                            iic=ir*pv.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*pv.ncol+ic];
                        //lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH,GlobalV::NLOCAL-i);
            Parallel_Reduce::reduce_complex_double_pool(lineS,GlobalV::NLOCAL-i);

            if (GlobalV::DRANK==0)
            {
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j-i],sizeof(std::complex<double>),1,g1);
                    fwrite(&lineS[j-i],sizeof(std::complex<double>),1,g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK==0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            for (int j=i; j<GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i*GlobalV::NLOCAL+j],sizeof(std::complex<double>),1,g1);
                fwrite(&S[i*GlobalV::NLOCAL+j],sizeof(std::complex<double>),1,g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } //end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK==0)
        {
            g1.open(ssh.str().c_str());
            g2.open(sss.str().c_str());
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir,ic;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            std::complex<double>* lineH = new std::complex<double>[GlobalV::NLOCAL-i];
            std::complex<double>* lineS = new std::complex<double>[GlobalV::NLOCAL-i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);

            ir = pv.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    ic = pv.trace_loc_col[j];
                    if (ic>=0)
                    {
                        int iic;
                        if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                        {
                            iic=ir+ic*pv.nrow;
                        }
                        else
                        {
                            iic=ir*pv.ncol+ic;
                        }
                        //lineH[j-i] = H[ir*pv.ncol+ic];
                        //lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j-i] = H[iic];
                        lineS[j-i] = S[iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH,GlobalV::NLOCAL-i);
            Parallel_Reduce::reduce_complex_double_pool(lineS,GlobalV::NLOCAL-i);

            if (GlobalV::DRANK==0)
            {
                for (int j=i; j<GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j-i];
                    g2 << " " << lineS[j-i];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        //if (GlobalV::DRANK==0);
        if (GlobalV::DRANK==0)           // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }

/*LiuXH add 2015-12-17,begin
        //int nprocs,myid;
        //MPI_Status status;
        //MPI_Comm_size(DIAG_HPSEPS_WORLD,&nprocs);
        //MPI_Comm_rank(DIAG_HPSEPS_WORLD,&myid);

        std::string H_fn;
        std::stringstream H_fn2;
        H_fn2<< "data-H-"  << GlobalV::DRANK ;
        H_fn=H_fn2.str();
        std::ofstream ofs_H;
        ofs_H.open(H_fn.c_str());
        ofs_H<<std::setprecision(8) << std::setw(12);

        std::string S_fn;
        std::stringstream S_fn2;
        S_fn2<< "data-S-"  << GlobalV::DRANK ;
        S_fn=S_fn2.str();
        std::ofstream ofs_S;
        ofs_S.open(S_fn.c_str());
        ofs_S<<std::setprecision(8) << std::setw(12);

        int irr,icc;
        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            irr = pv.trace_loc_row[i];
            if (irr>=0)
            {
                // data collection
                for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                        icc = pv.trace_loc_col[j];
                        if (icc>=0)
                        {
                                //if(abs(H[irr*pv.ncol+icc]) < 1.0e-10) H[irr*pv.ncol+icc] = 0.0;
                                //if(abs(S[irr*pv.ncol+icc]) < 1.0e-10) S[irr*pv.ncol+icc] = 0.0;
                                ofs_H << " " << H[irr*pv.ncol+icc];
                                ofs_S << " " << S[irr*pv.ncol+icc];
                        }
                }
                ofs_H << std::endl;
                ofs_S << std::endl;
             }
         }
//LiuXH add 2015-12-17,end*/
#else
        std::ofstream g1(ssh.str().c_str());
        std::ofstream g2(sss.str().c_str());

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i=0; i<GlobalV::NLOCAL; i++)
        {
            for (int j=i; j<GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i*GlobalV::NLOCAL+j];
                g2 << " " << S[i*GlobalV::NLOCAL+j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("HS_Matrix","save_HS_bit");
    return;
}

//void HS_Matrix::save_HSR_tr(const int Rx, const int Ry, const int Rz, const double *H, const double *S)
void HS_Matrix::save_HSR_tr(const int current_spin, LCAO_Matrix &lm)
//void HS_Matrix::save_HSR_tr(void)
{
    ModuleBase::TITLE("HS_Matrix","save_HSR_tr");
    ModuleBase::timer::tick("HS_Matrix","save_HSR_tr");

    std::stringstream ssh;
    std::stringstream sss;

    ssh << GlobalV::global_out_dir << "data-HR-tr_SPIN"<<current_spin;
    sss << GlobalV::global_out_dir << "data-SR-tr_SPIN"<<current_spin;
    //ssh << GlobalV::global_out_dir << "data-HR-tr_SPIN";
    //sss << GlobalV::global_out_dir << "data-SR-tr_SPIN";

#ifdef __MPI
    std::ofstream g1;
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        g1.open(ssh.str().c_str());
        g2.open(sss.str().c_str());
        g1 << "Matrix Dimension of H(R): "<<GlobalV::NLOCAL<<std::endl;
        g2 << "Matrix Dimension of S(R): "<<GlobalV::NLOCAL<<std::endl;
    }

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

//std::cout<<"R_x: "<<R_x<<std::endl;
//std::cout<<"R_y: "<<R_y<<std::endl;
//std::cout<<"R_z: "<<R_z<<std::endl;

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

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
//std::cout<<"dRx: "<<dRx<<std::endl;
//std::cout<<"dRy: "<<dRy<<std::endl;
//std::cout<<"dRz: "<<dRz<<std::endl;
                int ir,ic;
                for(int i=0; i<GlobalV::NLOCAL; i++)
                {
                    //double* lineH = new double[GlobalV::NLOCAL-i];
                    //double* lineS = new double[GlobalV::NLOCAL-i];
                    double* lineH = nullptr;
                    double* lineS = nullptr;
                    std::complex<double>* lineH_soc = nullptr;
                    std::complex<double>* lineS_soc = nullptr;
                    if(GlobalV::NSPIN!=4)
                    {
                        lineH = new double[GlobalV::NLOCAL];
                        lineS = new double[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
                        ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);
                    }
                    else
                    {
                        lineH_soc = new std::complex<double>[GlobalV::NLOCAL];
                        lineS_soc = new std::complex<double>[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(lineH_soc, GlobalV::NLOCAL);
                        ModuleBase::GlobalFunc::ZEROS(lineS_soc, GlobalV::NLOCAL);
                    }
                    //ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
                    //ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);
                    //ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
                    //ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

                    ir = lm.ParaV->trace_loc_row[i];
                    if(ir>=0)
                    {
                        //for(int j=i; j<GlobalV::NLOCAL; j++)
                        for(int j=0; j<GlobalV::NLOCAL; j++)
                        {
                            ic = lm.ParaV->trace_loc_col[j];
                            if(ic>=0)
                            {
                                //lineH[j-i] = H[ir*lm.ParaV->ncol+ic];
                                //lineS[j-i] = S[ir*lm.ParaV->ncol+ic];
                                int iic;
                                if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
                                {
                                    iic=ir+ic*lm.ParaV->nrow;
                                }
                                else
                                {
                                    iic=ir*lm.ParaV->ncol+ic;
                                }
                                if(GlobalV::NSPIN!=4)
                                {
                                    lineH[j] = lm.HR_tr[ix][iy][iz][iic];
                                    lineS[j] = lm.SlocR_tr[ix][iy][iz][iic];
                                }
                                else
                                {
                                    lineH_soc[j] = lm.HR_tr_soc[ix][iy][iz][iic];
                                    lineS_soc[j] = lm.SlocR_tr_soc[ix][iy][iz][iic];
                                }
                            }
                        }
                    }
                    else
                    {
                        //do nothing
                    }

                    //Parallel_Reduce::reduce_double_all(lineH,GlobalV::NLOCAL-i);
                    //Parallel_Reduce::reduce_double_all(lineS,GlobalV::NLOCAL-i);
                    if(GlobalV::NSPIN!=4)
                    {
                        Parallel_Reduce::reduce_double_all(lineH,GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_double_all(lineS,GlobalV::NLOCAL);
                    }
                    else
                    {
                        Parallel_Reduce::reduce_complex_double_all(lineH_soc,GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_complex_double_all(lineS_soc,GlobalV::NLOCAL);
                    }

                    if(GlobalV::DRANK==0)
                    {
                        //for(int j=i; j<GlobalV::NLOCAL; j++)
                        for(int j=0; j<GlobalV::NLOCAL; j++)
                        {
                            if(i==0 && j==0)
                            {
                                g1 << dRx << " " << dRy << " " << dRz  << "    //R std::vector(R2 - R1,unit: lattice vector)" <<std::endl;
                                g2 << dRx << " " << dRy << " " << dRz  << "    //R std::vector(R2 - R1,unit: lattice vector)" <<std::endl;
                            }
                            //g1 << " " << lineH[j-i];
                            //g2 << " " << lineS[j-i];
                            if(GlobalV::NSPIN!=4)
                            {
                                if(abs(lineH[j]) < 1.0e-12) lineH[j]=0.0;
                                if(abs(lineS[j]) < 1.0e-12) lineS[j]=0.0;
                                g1 << " " << lineH[j];
                                g2 << " " << lineS[j];
                            }
                            else
                            {
                if(abs(lineH_soc[j].real()) < 1.0e-12) lineH_soc[j]= std::complex<double> (0.0, lineH_soc[j].imag());
                                if(abs(lineH_soc[j].imag()) < 1.0e-12) lineH_soc[j]= std::complex<double> (lineH_soc[j].real(), 0.0);
                                if(abs(lineS_soc[j].real()) < 1.0e-12) lineS_soc[j]= std::complex<double> (0.0, lineS_soc[j].imag());
                                if(abs(lineS_soc[j].imag()) < 1.0e-12) lineS_soc[j]= std::complex<double> (lineS_soc[j].real() , 0.0);
                                g1 << " " << lineH_soc[j];
                                g2 << " " << lineS_soc[j];
                            }
                        }
                        g1 << std::endl;
                        g2 << std::endl;
                    }
                    if(GlobalV::NSPIN!=4)
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
                if(GlobalV::DRANK==0);
                {
                    g1.close();
                    g2.close();
                }
*/
            }
        }
    }
    //if(GlobalV::DRANK==0);
    if(GlobalV::DRANK==0)                // Peize Lin delete ; at 2020.01.31
    {
        g1.close();
        g2.close();
    }

#else
    std::ofstream g1(ssh.str().c_str());
    std::ofstream g2(sss.str().c_str());

    g1 << GlobalV::NLOCAL;
    g2 << GlobalV::NLOCAL;

    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        for (int j=i; j<GlobalV::NLOCAL; j++)
        {
            // not correct for serial version; need change 
            //g1 << " " << H[i*GlobalV::NLOCAL+j];
            //g2 << " " << S[i*GlobalV::NLOCAL+j];
        }
    }
    g1.close();
    g2.close();
#endif

    ModuleBase::timer::tick("HS_Matrix","save_HSR_tr");
    return;
}

void HS_Matrix::save_HSR_sparse(
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary,  
    const std::string &SR_filename, 
    const std::string &HR_filename_up, 
    const std::string &HR_filename_down = ""
)
{
    ModuleBase::TITLE("HS_Matrix","save_HSR_sparse");
    ModuleBase::timer::tick("HS_Matrix","save_HSR_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &HR_sparse_ptr = lm.HR_sparse;
    auto &SR_sparse_ptr = lm.SR_sparse;
    auto &HR_soc_sparse_ptr = lm.HR_soc_sparse;
    auto &SR_soc_sparse_ptr = lm.SR_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *H_nonzero_num[2] = {nullptr, nullptr};
    int *S_nonzero_num = nullptr;

    S_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(S_nonzero_num, total_R_num);

    int spin_loop = 1;
    if (GlobalV::NSPIN == 2)
    {
        spin_loop = 2;
    }

    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        H_nonzero_num[ispin] = new int[total_R_num];
        ModuleBase::GlobalFunc::ZEROS(H_nonzero_num[ispin], total_R_num);
    }

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                auto iter = HR_sparse_ptr[ispin].find(R_coor);
                if (iter != HR_sparse_ptr[ispin].end())
                {
                    for (auto &row_loop : iter->second)
                    {
                        H_nonzero_num[ispin][count] += row_loop.second.size();
                    }
                }
            }

            auto iter = SR_sparse_ptr.find(R_coor);
            if (iter != SR_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }
        else
        {
            auto iter = HR_soc_sparse_ptr.find(R_coor);
            if (iter != HR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    H_nonzero_num[0][count] += row_loop.second.size();
                }
            }

            iter = SR_soc_sparse_ptr.find(R_coor);
            if (iter != SR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(S_nonzero_num, total_R_num);
    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        Parallel_Reduce::reduce_int_all(H_nonzero_num[ispin], total_R_num);
    }

    if (GlobalV::NSPIN == 2)
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (H_nonzero_num[0][index] != 0 || H_nonzero_num[1][index] != 0 || S_nonzero_num[index] != 0)
            {
                output_R_number++;
            }
        }
    }
    else
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (H_nonzero_num[0][index] != 0 || S_nonzero_num[index] != 0)
            {
                output_R_number++;
            }
        }
    }

    std::stringstream ssh[2];
    std::stringstream sss;
    ssh[0] << GlobalV::global_out_dir << HR_filename_up;
    ssh[1] << GlobalV::global_out_dir << HR_filename_down;
    sss << GlobalV::global_out_dir << SR_filename;
    std::ofstream g1[2];
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                g1[ispin].open(ssh[ispin].str().c_str(), ios::binary);
                g1[ispin].write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
                g1[ispin].write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
            }

            g2.open(sss.str().c_str(), ios::binary);
            g2.write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
            g2.write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
        }
        else
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                g1[ispin].open(ssh[ispin].str().c_str());
                g1[ispin] << "Matrix Dimension of H(R): " << GlobalV::NLOCAL <<std::endl;
                g1[ispin] << "Matrix number of H(R): " << output_R_number << std::endl;
            }

            g2.open(sss.str().c_str());
            g2 << "Matrix Dimension of S(R): " << GlobalV::NLOCAL <<std::endl;
            g2 << "Matrix number of S(R): " << output_R_number << std::endl;
        }
    }

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (GlobalV::NSPIN == 2)
        {
            if (H_nonzero_num[0][count] == 0 && H_nonzero_num[1][count] == 0 && S_nonzero_num[count] == 0)
            {
                count++;
                continue;
            }
        }
        else
        {
            if (H_nonzero_num[0][count] == 0 && S_nonzero_num[count] == 0)
            {
                count++;
                continue;
            }
        }

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1[ispin].write(reinterpret_cast<char *>(&dRx), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&dRy), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&dRz), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&H_nonzero_num[ispin][count]), sizeof(int));
                }

                g2.write(reinterpret_cast<char *>(&dRx), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRy), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRz), sizeof(int));
                g2.write(reinterpret_cast<char *>(&S_nonzero_num[count]), sizeof(int));
            }
            else
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1[ispin] << dRx << " " << dRy << " " << dRz << " " << H_nonzero_num[ispin][count] << std::endl;
                }
                g2 << dRx << " " << dRy << " " << dRz << " " << S_nonzero_num[count] << std::endl;
            }
        }

        for (int ispin = 0; ispin < spin_loop; ++ispin)
        {
            if (H_nonzero_num[ispin][count] == 0)
            {
                // if (GlobalV::DRANK == 0)
                // {
                //     if (!binary)
                //     {
                //         g1[ispin] << std::endl;
                //         g1[ispin] << std::endl;
                //         for (int index = 0; index < GlobalV::NLOCAL+1; ++index)
                //         {
                //             g1[ispin] << 0 << " ";
                //         }
                //         g1[ispin] << std::endl;
                //     }
                // }
            }
            else
            {
                if (GlobalV::NSPIN != 4)
                {
                    output_single_R(g1[ispin], HR_sparse_ptr[ispin][R_coor], sparse_threshold, binary, *lm.ParaV);
                }
                else
                {
                    output_soc_single_R(g1[ispin], HR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
                }
            }
        }

        if (S_nonzero_num[count] == 0)
        {
            // if (!binary)
            // {
            //     if (GlobalV::DRANK == 0)
            //     {
            //         g2 << std::endl;
            //         g2 << std::endl;
            //         for (int index = 0; index < GlobalV::NLOCAL+1; ++index)
            //         {
            //             g2 << 0 << " ";
            //         }
            //         g2 << std::endl;
            //     }
            // }
        }
        else
        {
            if (GlobalV::NSPIN != 4)
            {
                output_single_R(g2, SR_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
            }
            else
            {
                output_soc_single_R(g2, SR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
            }
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        for (int ispin = 0; ispin < spin_loop; ++ispin) g1[ispin].close();
        g2.close();
    }
    
    for (int ispin = 0; ispin < spin_loop; ++ispin) 
    {
        delete[] H_nonzero_num[ispin];
        H_nonzero_num[ispin] = nullptr;
    }
    delete[] S_nonzero_num;
    S_nonzero_num = nullptr;

    ModuleBase::timer::tick("HS_Matrix","save_HSR_sparse");
    return;
}

void HS_Matrix::save_SR_sparse(
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary,  
    const std::string &SR_filename
)
{
    ModuleBase::TITLE("HS_Matrix","save_SR_sparse");
    ModuleBase::timer::tick("HS_Matrix","save_SR_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &SR_sparse_ptr = lm.SR_sparse;
    auto &SR_soc_sparse_ptr = lm.SR_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *S_nonzero_num = nullptr;

    S_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(S_nonzero_num, total_R_num);

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            auto iter = SR_sparse_ptr.find(R_coor);
            if (iter != SR_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }
        else
        {
            auto iter = SR_soc_sparse_ptr.find(R_coor);
            if (iter != SR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(S_nonzero_num, total_R_num);

    for (int index = 0; index < total_R_num; ++index)
    {
        if (S_nonzero_num[index] != 0)
        {
            output_R_number++;
        }
    }

    std::stringstream sss;
    sss << SR_filename;
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            g2.open(sss.str().c_str(), ios::binary);
            g2.write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
            g2.write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
        }
        else
        {
            g2.open(sss.str().c_str());
            g2 << "Matrix Dimension of S(R): " << GlobalV::NLOCAL <<std::endl;
            g2 << "Matrix number of S(R): " << output_R_number << std::endl;
        }
    }

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (S_nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                g2.write(reinterpret_cast<char *>(&dRx), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRy), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRz), sizeof(int));
                g2.write(reinterpret_cast<char *>(&S_nonzero_num[count]), sizeof(int));
            }
            else
            {
                g2 << dRx << " " << dRy << " " << dRz << " " << S_nonzero_num[count] << std::endl;
            }
        }

        if (GlobalV::NSPIN != 4)
        {
            output_single_R(g2, SR_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }
        else
        {
            output_soc_single_R(g2, SR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        g2.close();
    }

    delete[] S_nonzero_num;
    S_nonzero_num = nullptr;

    ModuleBase::timer::tick("HS_Matrix","save_SR_sparse");
    return;
}

void HS_Matrix::output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv)
{
    double *line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(GlobalV::NLOCAL + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << GlobalV::global_out_dir << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    line = new double[GlobalV::NLOCAL];
    for(int row = 0; row < GlobalV::NLOCAL; ++row)
    {
        // line = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);

        if(pv.trace_loc_row[row] >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (auto &value : iter->second)
                {
                    line[value.first] = value.second;
                }
            }
        }

        Parallel_Reduce::reduce_double_all(line, GlobalV::NLOCAL);

        if(GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < GlobalV::NLOCAL; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    if (binary)
                    {
                        ofs.write(reinterpret_cast<char *>(&line[col]), sizeof(double));
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        ofs << " " << fixed << scientific << std::setprecision(8) << line[col];
                        ofs_tem1 << " " << col;
                    }

                    nonzeros_count++;

                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }

        // delete[] line;
        // line = nullptr;

    }

    delete[] line;
    line = nullptr;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
            }
        }
        else
        {
            ofs << std::endl;
            ofs_tem1 << std::endl;
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str());
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }

        std::remove(tem1.str().c_str());

    }

}

void HS_Matrix::output_soc_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, std::complex<double>>> &XR, const double &sparse_threshold, const bool &binary, const Parallel_Orbitals &pv)
{
    std::complex<double> *line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(GlobalV::NLOCAL + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << GlobalV::global_out_dir << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
    }

    line = new std::complex<double>[GlobalV::NLOCAL];
    for(int row = 0; row < GlobalV::NLOCAL; ++row)
    {
        // line = new std::complex<double>[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);

        if(pv.trace_loc_row[row] >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (auto &value : iter->second)
                {
                    line[value.first] = value.second;
                }
            }
        }

        Parallel_Reduce::reduce_complex_double_all(line, GlobalV::NLOCAL);

        if (GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < GlobalV::NLOCAL; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    if (binary)
                    {
                        ofs.write(reinterpret_cast<char *>(&line[col]), sizeof(std::complex<double>));
                        ofs_tem1.write(reinterpret_cast<char *>(&col), sizeof(int));
                    }
                    else
                    {
                        ofs << " (" << fixed << scientific << std::setprecision(8) << line[col].real() << "," 
                                    << fixed << scientific << std::setprecision(8) << line[col].imag() << ")";
                        ofs_tem1 << " " << col;
                    }

                    nonzeros_count++;

                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }

        // delete[] line;
        // line = nullptr;

    }

    delete[] line;
    line = nullptr;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str(), ios::binary);
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs.write(reinterpret_cast<char *>(&i), sizeof(int));
            }
        }
        else
        {
            ofs << std::endl;
            ofs_tem1 << std::endl;
            ofs_tem1.close();
            ifs_tem1.open(tem1.str().c_str());
            ofs << ifs_tem1.rdbuf();
            ifs_tem1.close();
            for (auto &i : indptr)
            {
                ofs << " " << i;
            }
            ofs << std::endl;
        }

        std::remove(tem1.str().c_str());
    }

} 