#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"

void LCAO_Deepks::save_h_mat(const double *h_mat_in,const int nloc)
{
    ModuleBase::TITLE("LCAO_Deepks", "save_h_mat");
    for(int i=0;i<nloc;i++)
    {
        h_mat[i]=h_mat_in[i];
    }  
}

void LCAO_Deepks::save_h_mat(const std::complex<double> *h_mat_in,const int nloc)//for multi-k,not used now
{

}

void LCAO_Deepks::collect_h_mat(const double *h_in,ModuleBase::matrix &h_out,const int nlocal)
{
    ModuleBase::TITLE("LCAO_Deepks", "collect_h_tot");
    //construct the total H matrix
#ifdef __MPI
    int ir,ic;
    for (int i=0; i<nlocal; i++)
    {
        double* lineH = new double[nlocal-i];
        ModuleBase::GlobalFunc::ZEROS(lineH, nlocal-i);

        ir = pv->global2local_row(i);
        if (ir>=0)
        {
            // data collection
            for (int j=i; j<nlocal; j++)
            {
                ic = pv->global2local_col(j);
                if (ic>=0)
                {
                    int iic;
                    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                    {
                        iic=ir+ic*pv->nrow;
                    }
                    else
                    {
                        iic=ir*pv->ncol+ic;
                    }
                    //lineH[j-i] = H[ir*pv.ncol+ic];
                    lineH[j-i] = h_in[iic];
                }
            }
        }
        else
        {
            //do nothing
        }

        Parallel_Reduce::reduce_double_all(lineH,nlocal-i);

        for (int j=i; j<nlocal; j++)
        {
            h_out(i,j)=lineH[j-i];
            h_out(j,i)=lineH[j-i];//H is a symmetric matrix
        }
        delete[] lineH;
    }
#else
    for (int i=0; i<nlocal; i++)
    {
        for (int j=i; j<nlocal; j++)
        {
            h_out(i,j)=lineH[i*nlocal+j];
            h_out(j,i)=lineH[i*nlocal+j];//H is a symmetric matrix
        }
    }
#endif
}


#endif