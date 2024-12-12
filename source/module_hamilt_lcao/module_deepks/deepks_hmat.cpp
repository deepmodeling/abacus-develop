#ifdef __DEEPKS

#include "module_parameter/parameter.h"
#include "LCAO_deepks.h"
#include "module_base/parallel_reduce.h"

template <typename TK, typename TH>
void DeePKS_domain::collect_h_mat(
        const Parallel_Orbitals &pv,
		const std::vector<std::vector<TK>>& h_in,
		std::vector<TH> &h_out,
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
            std::vector<TK> lineH(nlocal-i, TK(0.0));

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

template <typename TH>
void DeePKS_domain::check_h_mat(
		const std::vector<TH> &H,
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

template void DeePKS_domain::collect_h_mat<double, ModuleBase::matrix>(
        const Parallel_Orbitals &pv,
        const std::vector<std::vector<double>>& h_in,
        std::vector<ModuleBase::matrix> &h_out,
        const int nlocal,
        const int nks);

template void DeePKS_domain::collect_h_mat<std::complex<double>, ModuleBase::ComplexMatrix>(
        const Parallel_Orbitals &pv,
        const std::vector<std::vector<std::complex<double>>>& h_in,
        std::vector<ModuleBase::ComplexMatrix> &h_out,
        const int nlocal,
        const int nks);

template void DeePKS_domain::check_h_mat<ModuleBase::matrix>(
        const std::vector<ModuleBase::matrix> &H,
        const std::string &h_file,
        const int nlocal,
        const int nks);

template void DeePKS_domain::check_h_mat<ModuleBase::ComplexMatrix>(
        const std::vector<ModuleBase::ComplexMatrix> &H,
        const std::string &h_file,
        const int nlocal,
        const int nks);

#endif
