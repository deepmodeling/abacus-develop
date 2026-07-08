#include "unitcell.h"
#include "source_base/parallel_common.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_xc/exx_info.h" // use GlobalC::exx_info

#include <string>
#include <vector>

namespace unitcell
{
#if defined(__MPI) && defined(__EXX)
    // Broadcast a vector<string> from rank 0 to all ranks.
    // Replaces the former cereal-based ModuleBase::bcast_data_cereal, which
    // was only ever used here to broadcast plain lists of ABFS file names and
    // pulled source_cell into a dependency on source_lcao/module_ri.
    static void bcast_string_vector(std::vector<std::string>& v)
    {
        int size = static_cast<int>(v.size());
        Parallel_Common::bcast_int(size);
        v.resize(size);
        for (int i = 0; i < size; ++i)
        {
            Parallel_Common::bcast_string(v[i]);
        }
    }
#endif

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
    
    void bcast_atoms_pseudo(Atom* atoms,
                                 const int ntype)
    {
    #ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = 0; i < ntype; i++) 
        {
            atoms[i].bcast_atom2();
        }
    #endif
    }

    void bcast_Lattice(Lattice& lat)
    {
    #ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        // distribute lattice parameters.
        ModuleBase::Matrix3& latvec = lat.latvec;
        ModuleBase::Matrix3& latvec_supercell = lat.latvec_supercell;
        Parallel_Common::bcast_string(lat.Coordinate);
        Parallel_Common::bcast_double(lat.lat0);
        Parallel_Common::bcast_double(lat.lat0_angstrom);
        Parallel_Common::bcast_double(lat.tpiba);
        Parallel_Common::bcast_double(lat.tpiba2);
        Parallel_Common::bcast_double(lat.omega);
        Parallel_Common::bcast_string(lat.latName);

        // distribute lattice vectors.
        Parallel_Common::bcast_double(latvec.e11);
        Parallel_Common::bcast_double(latvec.e12);
        Parallel_Common::bcast_double(latvec.e13);
        Parallel_Common::bcast_double(latvec.e21);
        Parallel_Common::bcast_double(latvec.e22);
        Parallel_Common::bcast_double(latvec.e23);
        Parallel_Common::bcast_double(latvec.e31);
        Parallel_Common::bcast_double(latvec.e32);
        Parallel_Common::bcast_double(latvec.e33);

         // distribute lattice vectors.
        for (int i = 0; i < 3; i++)
        {
            Parallel_Common::bcast_double(lat.a1[i]);
            Parallel_Common::bcast_double(lat.a2[i]);
            Parallel_Common::bcast_double(lat.a3[i]);
            Parallel_Common::bcast_double(lat.latcenter[i]);
            Parallel_Common::bcast_int(lat.lc[i]);
        }

        // distribute superlattice vectors.
        Parallel_Common::bcast_double(latvec_supercell.e11);
        Parallel_Common::bcast_double(latvec_supercell.e12);
        Parallel_Common::bcast_double(latvec_supercell.e13);
        Parallel_Common::bcast_double(latvec_supercell.e21);
        Parallel_Common::bcast_double(latvec_supercell.e22);
        Parallel_Common::bcast_double(latvec_supercell.e23);
        Parallel_Common::bcast_double(latvec_supercell.e31);
        Parallel_Common::bcast_double(latvec_supercell.e32);
        Parallel_Common::bcast_double(latvec_supercell.e33);

        // distribute Change the lattice vectors or not
    #endif
    }
    
    void bcast_magnetism(Magnetism& magnet, const int ntype)
    {
    #ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
        Parallel_Common::bcast_double(magnet.start_mag, ntype);
        if (PARAM.inp.nspin == 4) 
        {
            Parallel_Common::bcast_double(magnet.ux_[0]);
            Parallel_Common::bcast_double(magnet.ux_[1]);
            Parallel_Common::bcast_double(magnet.ux_[2]);
        }
    #endif
    }

    void bcast_unitcell(UnitCell& ucell)
    {
    #ifdef __MPI
        const int ntype = ucell.ntype;
        Parallel_Common::bcast_int(ucell.nat);

        bcast_Lattice(ucell.lat);
        bcast_magnetism(ucell.magnet,ntype);
        bcast_atoms_tau(ucell.atoms,ntype);

        for (int i = 0; i < ntype; i++)
        {
            Parallel_Common::bcast_string(ucell.orbital_fn[i]);
        }

        #ifdef __EXX
        bcast_string_vector(GlobalC::exx_info.info_ri.files_abfs);
        bcast_string_vector(GlobalC::exx_info.info_opt_abfs.files_abfs);
        bcast_string_vector(GlobalC::exx_info.info_opt_abfs.files_jles);
        #endif
        return;
    #endif
    }
}
