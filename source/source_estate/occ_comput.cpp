#include "occ_comput.h"

#include "source_base/tool_quit.h"

namespace elecstate
{

void occ_from_proj(
    const std::complex<double>* proj,
    const double* wg_ik,
    const int nbands,
    const int npol,
    const int nkb,
    const int nspin,
    const int isk,
    const int* nh_iat,
    const int nat,
    std::complex<double>* occ_block)
{
    if (proj == nullptr)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "proj is nullptr");
    }
    if (wg_ik == nullptr)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "wg_ik is nullptr");
    }
    if (nh_iat == nullptr)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "nh_iat is nullptr");
    }
    if (occ_block == nullptr)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "occ_block is nullptr");
    }
    if (nbands <= 0 || nkb <= 0 || nat <= 0)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "nbands, nkb and nat must be positive");
    }
    if (npol != 1 && npol != 2)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "npol must be 1 or 2");
    }
    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "nspin must be 1, 2 or 4");
    }
    if (nspin == 2 && isk != 0 && isk != 1)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "isk must be 0 or 1 when nspin=2");
    }
    if (nspin == 4 && npol != 2)
    {
        ModuleBase::WARNING_QUIT("occ_from_proj", "nspin=4 requires npol=2");
    }

    // rho^{ss'}_{iprj} = sum_i w_{k,i} * conj(proj^s_{i,iprj}) * proj^{s'}_{i,iprj}
    // iprj is the global projector index (begin_iprj + projector within atom)
    for (int ib = 0; ib < nbands; ib++)
    {
        const double weight = wg_ik[ib];
        int begin_iprj = 0;
        for (int iat = 0; iat < nat; iat++)
        {
            const int nprj = nh_iat[iat];
            for (int iprj = 0; iprj < nprj; iprj++)
            {
                const int occ_index = (begin_iprj + iprj) * 4;
                if (npol == 1)
                {
                    const int index = ib * nkb + begin_iprj + iprj;
                    const double occ = weight * (std::conj(proj[index]) * proj[index]).real();
                    if (nspin == 2 && isk == 1)
                    {
                        occ_block[occ_index + 3] += occ;
                    }
                    else if (nspin == 1)
                    {
                        // split evenly so the magnetization readout is zero
                        occ_block[occ_index] += 0.5 * occ;
                        occ_block[occ_index + 3] += 0.5 * occ;
                    }
                    else
                    {
                        occ_block[occ_index] += occ;
                    }
                }
                else
                {
                    // spinor components are offset by nkb in the proj layout
                    const int index = ib * 2 * nkb + begin_iprj + iprj;
                    occ_block[occ_index] += weight * std::conj(proj[index]) * proj[index];
                    occ_block[occ_index + 1] += weight * std::conj(proj[index]) * proj[index + nkb];
                    occ_block[occ_index + 2] += weight * std::conj(proj[index + nkb]) * proj[index];
                    occ_block[occ_index + 3] += weight * std::conj(proj[index + nkb]) * proj[index + nkb];
                }
            }
            begin_iprj += nprj;
        }
    }
}

} // namespace elecstate
