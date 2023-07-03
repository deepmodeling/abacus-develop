#include "paw_cell.h"
#include "module_base/tool_title.h"

void Paw_Cell::init_paw_cell(
    const double ecutwfc_in, const double cell_factor_in,
    const int nat_in, const int ntyp_in,
    const int * atom_type_in, const double ** atom_coord_in,
    const std::vector<std::string> filename_list_in,
    const int nx_in, const int ny_in, const int nz_in,
    const double * eigts1_in, const double * eigts2_in, const double * eigts3_in)
{
    ModuleBase::TITLE("Paw_Element","init_paw_cell");

    this -> nat = nat_in;
    this -> ntyp = ntyp_in;
    this -> nx = nx_in;
    this -> ny = ny_in;
    this -> nz = nz_in;

    atom_coord.resize(nat);
    atom_type.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        atom_coord[iat].resize(3);
        for(int i = 0; i < 3; i ++)
        {
            atom_coord[iat][i] = atom_coord_in[iat][i];
        }

        atom_type[iat] = atom_type_in[iat];
    }

    paw_element_list.resize(ntyp);
    assert(filename_list_in.size() == ntyp);
    for(int ityp = 0; ityp < ntyp; ityp ++)
    {
        paw_element_list[ityp].init_paw_element(ecutwfc_in, cell_factor_in);
        paw_element_list[ityp].read_paw_xml(filename_list_in[ityp]);
    }

    eigts1.resize(nat);
    eigts2.resize(nat);
    eigts3.resize(nat);

    for(int iat = 0; iat < nat; iat ++)
    {
        eigts1[iat].resize(2*nx+1);
        eigts2[iat].resize(2*ny+1);
        eigts3[iat].resize(2*nz+1);

        for(int i = 0; i < 2*nx+1; i ++)
        {
            eigts1[iat][i] = eigts1_in[ iat * (2*nx+1) + i];
        }

        for(int i = 0; i < 2*ny+1; i ++)
        {
            eigts2[iat][i] = eigts2_in[ iat * (2*ny+1) + i];
        }

        for(int i = 0; i < 2*nz+1; i ++)
        {
            eigts3[iat][i] = eigts3_in[ iat * (2*nz+1) + i];
        }
    }
}

// exp(-i(k+G)R_I) = exp(-ikR_I) exp(-iG_xR_Ix) exp(-iG_yR_Iy) exp(-iG_zR_Iz)

void Paw_Cell::set_paw_k(
    const int npw, const double * kpt,
    const int * ig_to_ix, const int * ig_to_iy, const int * ig_to_iz)
{
    ModuleBase::TITLE("Paw_Element","set_paw_k");

    const double pi = 3.141592653589793238462643383279502884197;
    const double twopi = 2.0 * pi;

    struc_fact.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        double arg = 0.0;
        for(int i = 0; i < 3; i ++)
        {
            arg += atom_coord[iat][i] * kpt[i];
        }
        arg *= twopi;
        const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));

        struc_fact[iat].resize(npw);

        for(int ipw = 0; ipw < npw; ipw ++)
        {
            int ix = ig_to_ix[ipw];
            int iy = ig_to_iy[ipw];
            int iz = ig_to_iz[ipw];

            struc_fact[iat][ipw] = kphase * eigts1[iat][ix] * eigts2[iat][iy] * eigts3[iat][iz];
        }
    }
}