#include "paw_atom.h"
#include "module_base/global_variable.h"

void Paw_Atom::init_paw_atom(const int nproj_in)
{

    nproj = nproj_in;

    ca.resize(nproj);

    rhoij.resize(GlobalV::NSPIN);
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        rhoij[is].resize(nproj*(nproj + 1) / 2);
    }

    rhoijp.resize(GlobalV::NSPIN * nproj*(nproj + 1) / 2);
    rhoijselect.resize(nproj*(nproj + 1) / 2);

    dij.resize(GlobalV::NSPIN);
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        dij[is].resize(nproj*nproj);
    }
    sij.resize(nproj*nproj);

    this -> reset_rhoij();
    this -> reset_dij();
    this -> reset_sij();
}

void Paw_Atom::set_ca(std::vector<std::complex<double>> & ca_in, const double weight_in)
{
    for(int i = 0; i < nproj; i ++)
    {
        ca[i] = ca_in[i];
    }

    weight = weight_in;
}

void Paw_Atom::reset_rhoij()
{
    nrhoijsel = 0;

    for(int i = 0; i < nproj*(nproj+1)/2; i ++)
    {
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            rhoij[is][i] = 0.0;
        }
        rhoijselect[i] = -1;
    }    

    for(int i = 0; i < GlobalV::NSPIN * nproj*(nproj + 1) / 2; i ++)
    {
        rhoijp[i] = 0.0;
    }
}

void Paw_Atom::accumulate_rhoij(const int current_spin)
{
    for(int iproj = 0; iproj < nproj; iproj ++)
    {
        int i0 = iproj * (iproj + 1) / 2;
        for(int jproj = 0; jproj < iproj+1; jproj ++)
        {
            std::complex<double> tmp = std::conj(ca[iproj]) * ca[jproj];
            rhoij[current_spin][i0 + jproj] += tmp.real() * weight;
        }
    }
}

void Paw_Atom::convert_rhoij()
{
    nrhoijsel = 0;
    for(int i = 0; i < rhoij[0].size(); i ++)
    {
        bool nonzero = false;
        for(int is = 0; is < GlobalV::NSPIN; is ++)
        {
            if(std::abs(rhoij[is][i]) > 1e-10)
            {
                nonzero = true;
                break;
            }
        }

        if(nonzero)
        {
            rhoijselect[nrhoijsel] = i+1; //index in fortran
            for(int is = 0; is < GlobalV::NSPIN; is ++)
            {
                rhoijp[nrhoijsel + is * rhoij[0].size()] = rhoij[is][i];
            }
            nrhoijsel ++;
        }
    }
}

void Paw_Atom::reset_dij()
{
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int i = 0; i < nproj*nproj; i ++)
        {
            dij[is][i] = 0.0;
        }
    }
}

void Paw_Atom::set_dij(double** dij_in)
{
    for(int is = 0; is < GlobalV::NSPIN; is ++)
    {
        for(int i = 0; i < nproj*nproj; i ++)
        {
            dij[is][i] = dij_in[is][i];
        }
    }
}

void Paw_Atom::reset_sij()
{
    for(int i = 0; i < nproj*nproj; i ++)
    {
        sij[i] = 0.0;
    }
}

void Paw_Atom::set_sij(const double* sij_in)
{
    for(int i = 0; i < nproj*nproj; i ++)
    {
        sij[i] = sij_in[i];
    }
}