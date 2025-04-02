#include "unitcell.h"   
#include "module_base/parallel_common.h"
#include "module_parameter/parameter.h"
#ifdef __EXX
#include "module_ri/serialization_cereal.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif
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
        Parallel_Common::bcast_double(magnet.start_magnetization, ntype);
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
        ModuleBase::bcast_data_cereal(GlobalC::exx_info.info_ri.files_abfs,
                                    MPI_COMM_WORLD,
                                    0);
        #endif
        return;
    #endif
    }
}