#include "module_parameter/parameter.h"

#ifdef __DEEPKS

#include "LCAO_deepks_io.h"
#include "module_base/tool_quit.h"
#include "npy.hpp"

#include <mpi.h>

template <typename TK>
void LCAO_deepks_io::print_dm(const int nks, const int nlocal, const int nrow, const std::vector<std::vector<TK>>& dm)
{
    std::stringstream ss;
    for (int ik = 0; ik < nks; ik++)
    {
        ss.str("");
        ss << "dm_" << ik;
        std::ofstream ofs(ss.str().c_str());
        ofs << std::setprecision(15);

        for (int mu = 0; mu < nlocal; mu++)
        {
            for (int nu = 0; nu < nlocal; nu++)
            {
                ofs << dm[ik][mu * nrow + nu] << " ";
            }
            ofs << std::endl;
        }
    }
}

void LCAO_deepks_io::load_npy_gedm(const int nat,
                                   const int des_per_atom,
                                   double** gedm,
                                   double& e_delta,
                                   const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "load_npy_gedm");

    if (rank == 0)
    {
        // load gedm.npy
        std::vector<double> npy_gedm;
        std::vector<unsigned long> dshape = {static_cast<unsigned long>(nat), static_cast<unsigned long>(des_per_atom)};

        std::string gedm_file = "gedm.npy";

        npy::LoadArrayFromNumpy(gedm_file, dshape, npy_gedm);

        for (int iat = 0; iat < nat; iat++)
        {
            for (int ides = 0; ides < des_per_atom; ides++)
            {
                gedm[iat][ides] = npy_gedm[iat * des_per_atom + ides] * 2.0; // Ha to Ry
            }
        }

        // load ec.npy
        std::vector<double> npy_ec;
        std::vector<unsigned long> eshape = {1ul};
        std::string ec_file = "ec.npy";
        npy::LoadArrayFromNumpy(ec_file, eshape, npy_ec);
        e_delta = npy_ec[0] * 2.0; // Ha to Ry
    }

#ifdef __MPI
    for (int iat = 0; iat < nat; iat++)
    {
        MPI_Bcast(gedm[iat], des_per_atom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&e_delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

// saves descriptor into dm_eig.npy
void LCAO_deepks_io::save_npy_d(const int nat,
                                const int des_per_atom,
                                const int inlmax,
                                const int* inl_l,
                                const bool deepks_equiv,
                                const std::vector<torch::Tensor>& descriptor,
                                const std::string& out_dir,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_d");

    if (rank != 0)
    {
        return;
    }

    // save descriptor in .npy format
    //  deepks_equiv was PARAM.inp.deepks_equiv
    if (!deepks_equiv)
    {
        std::vector<double> npy_des;
        for (int inl = 0; inl < inlmax; ++inl)
        {
            auto accessor = descriptor[inl].accessor<double, 1>();
            int nm = 2 * inl_l[inl] + 1;
            for (int im = 0; im < nm; im++)
            {
                npy_des.push_back(accessor[im]);
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(des_per_atom)};
        if (rank == 0)
        {
            std::string file_dm_eig = out_dir + "deepks_dm_eig.npy";
            // std::string file_dm_eig = "dm_eig.npy";
            npy::SaveArrayAsNumpy(file_dm_eig, false, 2, dshape, npy_des);
        }
    }
    else
    {
        // a rather unnecessary way of writing this, but I'll do it for now
        std::vector<double> npy_des;
        for (int iat = 0; iat < nat; iat++)
        {
            auto accessor = descriptor[iat].accessor<double, 1>();
            for (int i = 0; i < des_per_atom; i++)
            {
                npy_des.push_back(accessor[i]);
            }
        }
        const long unsigned dshape[] = {static_cast<unsigned long>(nat), static_cast<unsigned long>(des_per_atom)};
        if (rank == 0)
        {
            std::string file_dm_eig = out_dir + "deepks_dm_eig.npy";
            // std::string file_dm_eig = "dm_eig.npy";
            npy::SaveArrayAsNumpy(file_dm_eig, false, 2, dshape, npy_des);
        }
    }
    return;
}

// saves energy in numpy format
void LCAO_deepks_io::save_npy_e(const double& e, const std::string& e_file, const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_e");
    if (rank != 0)
    {
        return;
    }

    // save energy in .npy format
    const long unsigned eshape[] = {1};
    npy::SaveArrayAsNumpy(e_file, false, 1, eshape, &e);
    return;
}

template <typename TK, typename TH>
void LCAO_deepks_io::save_npy_h(const std::vector<TH>& hamilt,
                                const std::string& h_file,
                                const int nlocal,
                                const int nks,
                                const int rank)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_npy_h");
    if (rank != 0)
    {
        return;
    }

    const long unsigned hshape[]
        = {static_cast<unsigned long>(nks), static_cast<unsigned long>(nlocal), static_cast<unsigned long>(nlocal)};

    std::vector<TK> npy_h;
    for (int k = 0; k < nks; k++)
    {
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                npy_h.push_back(hamilt[k](i, j));
            }
        }
    }

    npy::SaveArrayAsNumpy(h_file, false, 3, hshape, npy_h);
    return;
}

void LCAO_deepks_io::save_matrix2npy(const std::string& file_name,
                                     const ModuleBase::matrix& matrix,
                                     const int rank,
                                     const double& scale,
                                     const char mode)
{
    ModuleBase::TITLE("LCAO_deepks_io", "save_matrix2npy");

    if (rank != 0)
    {
        return;
    }
    const int nr = matrix.nr;
    const int nc = matrix.nc;
    int size = 0;
    std::vector<long unsigned> shape;

    if (mode == 'U' || mode == 'L') // upper or lower triangular
    {
        assert(nr == nc);
        size = nr * (nr + 1) / 2;
        shape.resize(1);
        shape[0] = size;
    }
    else if (mode == 'N') // normal
    {
        size = nr * nc;
        shape.resize(2);
        shape[0] = nr;
        shape[1] = nc;
    }
    else
    {
        ModuleBase::WARNING_QUIT("save_matrix2npy", "Invalid mode! Support only 'U', 'L', 'N'.");
    }

    std::vector<double> scaled_data(size);
    if (mode == 'U') // upper triangular
    {
        int index = 0;
        for (int i = 0; i < nr; ++i)
        {
            for (int j = i; j < nc; ++j)
            {
                scaled_data[index] = matrix(i, j) * scale;
                index++;
            }
        }
    }
    else if (mode == 'L') // lower triangular
    {
        int index = 0;
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                scaled_data[index] = matrix(i, j) * scale;
                index++;
            }
        }
    }
    else // normal
    {
        for (int i = 0; i < nr; ++i)
        {
            for (int j = 0; j < nc; ++j)
            {
                scaled_data[i * nc + j] = matrix(i, j) * scale;
            }
        }
    }

    npy::SaveArrayAsNumpy(file_name, false, shape.size(), shape.data(), scaled_data);
    return;
}

template <typename T>
void LCAO_deepks_io::save_tensor2npy(const std::string& file_name, const torch::Tensor& tensor, const int rank)
{
    if (rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("LCAO_deepks_io", "save_tensor2npy");
    const int dim = tensor.dim();
    std::vector<long unsigned> shape(dim);
    for (int i = 0; i < dim; i++)
    {
        shape[i] = tensor.size(i);
    }

    std::vector<T> data(tensor.numel());

    if constexpr (std::is_same<T, double>::value)
    {
        std::memcpy(data.data(), tensor.data_ptr<double>(), tensor.numel() * sizeof(double));
    }
    else
    {
        auto tensor_data = tensor.data_ptr<c10::complex<double>>();
        for (size_t i = 0; i < tensor.numel(); ++i)
        {
            data[i] = std::complex<double>(tensor_data[i].real(), tensor_data[i].imag());
        }
    }

    npy::SaveArrayAsNumpy(file_name, false, shape.size(), shape.data(), data);
}

template void LCAO_deepks_io::print_dm<double>(const int nks,
                                               const int nlocal,
                                               const int nrow,
                                               const std::vector<std::vector<double>>& dm);

template void LCAO_deepks_io::print_dm<std::complex<double>>(const int nks,
                                                             const int nlocal,
                                                             const int nrow,
                                                             const std::vector<std::vector<std::complex<double>>>& dm);

template void LCAO_deepks_io::save_npy_h<double>(const std::vector<ModuleBase::matrix>& hamilt,
                                                 const std::string& h_file,
                                                 const int nlocal,
                                                 const int nks,
                                                 const int rank);

template void LCAO_deepks_io::save_npy_h<std::complex<double>>(const std::vector<ModuleBase::ComplexMatrix>& hamilt,
                                                               const std::string& h_file,
                                                               const int nlocal,
                                                               const int nks,
                                                               const int rank);

template void LCAO_deepks_io::save_tensor2npy<double>(const std::string& file_name,
                                                      const torch::Tensor& tensor,
                                                      const int rank);

template void LCAO_deepks_io::save_tensor2npy<std::complex<double>>(const std::string& file_name,
                                                                    const torch::Tensor& tensor,
                                                                    const int rank);
#endif
