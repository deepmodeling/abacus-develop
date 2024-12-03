#ifdef __DEEPKS

#include "module_parameter/parameter.h"
#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"


void DeePKS_domain::save_h_mat(
		const double *h_mat_in,
		const int nloc,
        const int ik)
{
	for(int i=0;i<nloc;i++)
	{
		GlobalC::ld.h_mat[i]=h_mat_in[i];
	}  
}

void DeePKS_domain::save_h_mat(
		const std::complex<double> *h_mat_in,
		const int nloc,
        const int ik)
{
    for(int i=0;i<nloc;i++)
    {
        GlobalC::ld.h_mat_k[ik][i]=h_mat_in[i];
    }  
}


void DeePKS_domain::collect_h_mat(
        const Parallel_Orbitals &pv,
		const std::vector<double>& h_in,
		ModuleBase::matrix &h_out,
		const int nlocal)
{
    ModuleBase::TITLE("DeePKS_domain", "collect_h_tot");

    //construct the total H matrix
#ifdef __MPI
    int ir=0;
    int ic=0;
    for (int i=0; i<nlocal; i++)
    {
        std::vector<double> lineH(nlocal-i,0.0);

        ir = pv.global2local_row(i);
        if (ir>=0)
        {
            // data collection
            for (int j=i; j<nlocal; j++)
            {
                ic = pv.global2local_col(j);
                if (ic>=0)
                {
                    int iic=0;
                    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                    {
                        iic=ir+ic*pv.nrow;
                    }
                    else
                    {
                        iic=ir*pv.ncol+ic;
                    }
                    lineH[j-i] = h_in[iic];
                }
            }
        }
        else
        {
            //do nothing
        }

        Parallel_Reduce::reduce_all(lineH.data(),nlocal-i);

        for (int j=i; j<nlocal; j++)
        {
            h_out(i,j)=lineH[j-i];
            h_out(j,i)=h_out(i,j);//H is a symmetric matrix
        }
    }
#else
    for (int i=0; i<nlocal; i++)
    {
        for (int j=i; j<nlocal; j++)
        {
            h_out(i,j)=h_in[i*nlocal+j];
            h_out(j,i)=h_out(i,j);//H is a symmetric matrix
        }
    }
#endif
}

// for multi-k
void DeePKS_domain::collect_h_mat(
        const Parallel_Orbitals &pv,
		const std::vector<std::vector<std::complex<double>>>& h_in,
		std::vector<ModuleBase::ComplexMatrix> &h_out,
		const int nlocal,
        const int nks)
{
    ModuleBase::TITLE("DeePKS_domain", "collect_h_tot");

    //construct the total H matrix
    for (int k=0; k<nks; k++) 
    {
#ifdef __MPI
        int ir=0;
        int ic=0;
        for (int i=0; i<nlocal; i++)
        {
            std::vector<std::complex<double>> lineH(nlocal-i, std::complex<double>(0.0, 0.0));

            ir = pv.global2local_row(i);
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<nlocal; j++)
                {
                    ic = pv.global2local_col(j);
                    if (ic>=0)
                    {
                        int iic=0;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
                        {
                            iic=ir+ic*pv.nrow;
                        }
                        else
                        {
                            iic=ir*pv.ncol+ic;
                        }
                        lineH[j-i] = h_in[k][iic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_all(lineH.data(),nlocal-i);

            for (int j=i; j<nlocal; j++)
            {
                h_out[k](i,j)=lineH[j-i];
                h_out[k](j,i)=h_out[k](i,j);//H is a symmetric matrix
            }
        }
#else
        for (int i=0; i<nlocal; i++)
        {
            for (int j=i; j<nlocal; j++)
            {
                h_out[k](i,j)=h_in[k][i*nlocal+j];
                h_out[k](j,i)=h_out[k](i,j);//H is a symmetric matrix
            }
        }
#endif
    }
}

//just for gamma-only now
void DeePKS_domain::check_h_mat(
		const ModuleBase::matrix &H,
		const std::string &h_file,
		const int nlocal)
{
    std::ofstream ofs(h_file.c_str());
    ofs << std::setprecision(10);
    for (int i=0; i<nlocal; i++)
    {
        for (int j=0; j<nlocal; j++)
        {
            ofs << H(i,j) << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}

// for multi-k
void DeePKS_domain::check_h_mat(
		const std::vector<ModuleBase::ComplexMatrix> &H,
		const std::string &h_file,
		const int nlocal,
        const int nks)
{
    std::ofstream ofs(h_file.c_str());
    ofs << std::setprecision(10);
    for (int k=0; k<nks; k++)
    {
        for (int i=0; i<nlocal; i++)
        {
            for (int j=0; j<nlocal; j++)
            {
                ofs << H[k](i,j) << " ";
            }
            ofs << std::endl;
        }
        ofs << std::endl;
    }
    ofs.close();
}

#endif
