#include "unitcell.h"   

namespace unitcell
{
    void bcast_atoms_tau(Atom* atoms,
                         const int ntype)
    {
    #ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < ntype; i++) {
            atoms[i].bcast_atom(); // bcast tau array
        }
    #endif
    }
}