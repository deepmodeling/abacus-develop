#include "../hcontainer.h"

void H_gamma_to_matrix(const hamilt::HContainer<double>& hR, 
    std::complex<double>* hk, 
    const int ncol)
{
    // loop_atom_pairs
    for (int j = 0; j < hR.size_atom_pair(); ++j)
    {
        // Hk += HR * e^ikR
        hR.get_atom_pair(j).add_to_matrix(hk, ncol, 1.0, 0);
    }
}