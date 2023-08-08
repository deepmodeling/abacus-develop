#include "paw_cell.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

// The subroutines here are used to gather information from the main ABACUS program
// 1. ecut, ecutpaw : kinetic energy cutoff of the planewave basis set
// there will be one coarse grid for density/potential, and a fine grid for PAW
// the unit is in Hartree
// 2. rprimd, gprimd : real and reciprocal space lattice vectors, respectively
// unit for rprimd is in Bohr, and for gprimd is in Bohr^-1
// 3. gmet : reciprocal space metric (bohr^-2)
// 4. ucvol : volume of unit cell (Bohr^3)
// 5. ngfft, ngfftdg : dimension of FFT grids of the corase and fine grids
// 6. natom, ntypat, typat: #. atoms, #. element types
// and typat records the type of each atom
// 7. xred : coordinate of each atom, in terms of rprimd (namely, direct coordinate)
// 8. filename_list : filename of the PAW xml files for each element

// Cutoff energies, sets ecut and ecutpaw
void Paw_Cell::set_libpaw_ecut(const double ecut_in, const double ecutpaw_in)
{
    ecut = ecut_in;
    ecutpaw = ecutpaw_in;
}

// Sets rprimd, gprimd, gmet and ucvol
// Only real space lattice vector needed, others are to be calculated
void Paw_Cell::set_libpaw_cell(const ModuleBase::Matrix3 latvec, const double lat0)
{
    rprimd.resize(9);
    gprimd.resize(9);
    gmet.resize(9);
}

// FFT grid information, sets ngfft and ngfftdg
void Paw_Cell::set_libpaw_fft(const int * ngfft_in, const int * ngfftdg_in)
{
    ngfft.resize(3);
    ngfftdg.resize(3);

    for(int i = 0; i < 3; i ++)
    {
        ngfft[i] = ngfft_in[i];
        ngfftdg[i] = ngfftdg_in[i];
    }
}

// Sets natom, ntypat, typat and xred
void Paw_Cell::set_libpaw_atom(const int natom_in, const int ntypat_in, const int * typat_in, const double * xred_in)
{
    natom = natom_in;
    ntypat = ntypat_in;

    typat.resize(natom);
    xred.resize(3*natom);
}

// Sets filename_list
void Paw_Cell::set_libpaw_files()
{

}