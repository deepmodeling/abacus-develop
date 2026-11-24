#ifndef DEEPKS_EINSUM_H
#define DEEPKS_EINSUM_H

#include <torch/torch.h>

namespace deepks
{

inline torch::Tensor einsum_kamn_avmn(const torch::Tensor& a, const torch::Tensor& b)
{
    const int k = a.size(0);
    const int na = a.size(1);
    const int v = b.size(1);
    auto out = torch::zeros({k, na, v}, a.options());
    for (int ik = 0; ik < k; ++ik)
    {
        for (int ia = 0; ia < na; ++ia)
        {
            for (int iv = 0; iv < v; ++iv)
            {
                out[ik][ia][iv] = (a[ik][ia] * b[ia][iv]).sum();
            }
        }
    }
    return out;
}

inline torch::Tensor einsum_kxyamn_avmn(const torch::Tensor& a, const torch::Tensor& b)
{
    const int k = a.size(0);
    const int nx = a.size(1);
    const int ny = a.size(2);
    const int na = a.size(3);
    const int v = b.size(1);
    auto out = torch::zeros({k, nx, ny, na, v}, a.options());
    for (int ik = 0; ik < k; ++ik)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ia = 0; ia < na; ++ia)
                {
                    for (int iv = 0; iv < v; ++iv)
                    {
                        out[ik][ix][iy][ia][iv] = (a[ik][ix][iy][ia] * b[ia][iv]).sum();
                    }
                }
            }
        }
    }
    return out;
}

inline torch::Tensor einsum_bxamn_avmn(const torch::Tensor& a, const torch::Tensor& b)
{
    const int nb = a.size(0);
    const int nx = a.size(1);
    const int na = a.size(2);
    const int v = b.size(1);
    auto out = torch::zeros({nb, nx, na, v}, a.options());
    for (int ib = 0; ib < nb; ++ib)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            for (int ia = 0; ia < na; ++ia)
            {
                for (int iv = 0; iv < v; ++iv)
                {
                    out[ib][ix][ia][iv] = (a[ib][ix][ia] * b[ia][iv]).sum();
                }
            }
        }
    }
    return out;
}

inline torch::Tensor einsum_bamn_avmn(const torch::Tensor& a, const torch::Tensor& b)
{
    const int nb = a.size(0);
    const int na = a.size(1);
    const int v = b.size(1);
    auto out = torch::zeros({nb, na, v}, a.options());
    for (int ib = 0; ib < nb; ++ib)
    {
        for (int ia = 0; ia < na; ++ia)
        {
            for (int iv = 0; iv < v; ++iv)
            {
                out[ib][ia][iv] = (a[ib][ia] * b[ia][iv]).sum();
            }
        }
    }
    return out;
}

} // namespace deepks

#endif
