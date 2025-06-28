#ifdef __MLALGO

#include "deepks_check.h"

template <typename T>
void DeePKS_domain::check_tensor(const torch::Tensor& tensor, const std::string& filename, const int rank)
{
    if (rank != 0)
    {
        return;
    }
    using T_tensor =
        typename std::conditional<std::is_same<T, std::complex<double>>::value, c10::complex<double>, T>::type;

    std::ofstream ofs(filename.c_str());
    ofs << std::setprecision(10);

    auto sizes = tensor.sizes();
    int ndim = sizes.size();
    auto data_ptr = tensor.data_ptr<T_tensor>();
    int64_t numel = tensor.numel();

    // stride for each dimension
    std::vector<int64_t> strides(ndim, 1);
    for (int i = ndim - 2; i >= 0; --i)
    {
        strides[i] = strides[i + 1] * sizes[i + 1];
    }

    for (int64_t idx = 0; idx < numel; ++idx)
    {
        // index to multi-dimensional indices
        std::vector<int64_t> indices(ndim);
        int64_t tmp = idx;
        for (int d = 0; d < ndim; ++d)
        {
            indices[d] = tmp / strides[d];
            tmp = tmp % strides[d];
        }

        T_tensor tmp_val = data_ptr[idx];
        T* tmp_ptr = reinterpret_cast<T*>(&tmp_val);
        ofs << *tmp_ptr;

        // print space or newline
        if (((idx + 1) % sizes[ndim - 1]) == 0)
        {
            ofs << std::endl;
        }
        else
        {
            ofs << " ";
        }
    }

    ofs.close();
}



template void DeePKS_domain::check_tensor<int>(const torch::Tensor& tensor, const std::string& filename, const int rank);
template void DeePKS_domain::check_tensor<double>(const torch::Tensor& tensor, const std::string& filename, const int rank);
template void DeePKS_domain::check_tensor<std::complex<double>>(const torch::Tensor& tensor, const std::string& filename, const int rank);

#endif
