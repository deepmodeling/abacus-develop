#ifndef PSEUDO_H
#define PSEUDO_H

#include "module_base/global_function.h"
#include "module_io/output.h"

//-----------------------------------------
// read in pseudopotentials
// mohan update 2021-05-01
//-----------------------------------------
class pseudo
{
  public:
    pseudo();
    ~pseudo();

    //<PP_HEADER>
    bool has_so = false;  // if .true. includes spin-orbit
    int nv = 0;           // UPF file version number
    std::string psd;      // Element label
    std::string pp_type;  // Pseudo type ( NC or US )
    bool tvanp = false;   // .true. if Ultrasoft
    bool nlcc = false;    // Non linear core corrections(bool)
    std::string xc_func;  // Exch-Corr type
    double zv = 0;           // z valence
    double etotps = 0.0;  // total energy
    double ecutwfc = 0.0; // suggested cut-off for wfc
    double ecutrho = 0.0; // suggested cut-off for rho
    int lmax = 0;         // maximum angular momentum component
    int mesh = 0;         // number of point in the radial mesh
    int nchi = 0;         // nwfc,number of wavefunctions
    int nbeta = 0;        // number of projectors
    int nqlc = 0;         // number of angular momenta in Q
    int kkbeta = 0;       // kkbeta, point where the beta are zero

    std::vector<std::string> els = {}; // els[nchi]
    std::vector<int> lchi = {};        // lchi[nchi]
    std::vector<double> oc = {};       // oc[nchi]

    std::vector<double> jjj = {};  // total angual momentum, jjj[nbeta]
    std::vector<double> jchi = {}; // jchi(nwfc), added by zhengdy-soc
    std::vector<int> nn = {};

    // Local pseudopotentials
    std::vector<double> vloc_at = {}; // [mesh], local potential( = pseudopot_upf.vloc )

    // <PP_MESH>
    std::vector<double> r = {};   // radial logaritmic mesh, r[0:mesh-1]
    std::vector<double> rab = {}; // derivative of the radial mesh, rab[0:mesh-1]

    //<PP_NLCC>
    std::vector<double> rho_atc = {}; // radial core charge density, rho_atc[0:mesh-1]

    //<PP_RHOATOM>
    std::vector<double> rho_at = {}; // radial atomic charge density, rho_at[0:mesh-1]

    // <PP_PSWFC>
    ModuleBase::matrix chi; // radial atomic orbitals, chi(nchi, mesh)

    // other
    int msh = 0;       // number of points up to rcut
    double rcut = 0.0; // cut-off radius

    // <PP_BETA>
    std::vector<int> lll = {}; // lll(nbeta), angular momentum of the beta function

    // <PP_DIJ>
    ModuleBase::matrix dion;  // dion(nbeta,nbeta)
    ModuleBase::matrix betar; // (nbeta, mesh), radial beta_{mu} functions

    // other
    int nh = 0; // number of beta functions per atomic type

    // uspp
    ModuleBase::realArray qfuncl; // qfuncl(2*lmax+1,nbeta*(nbeta+1)/2,mesh) Q_{mu,nu}(|r|) function for |r|> r_L
    ModuleBase::matrix qqq;       // qqq(nbeta,nbeta) q_{mu,nu}

    void print_pseudo_h(std::ofstream& ofs);
    void print_pseudo_atom(std::ofstream& ofs);
    void print_pseudo_vl(std::ofstream& ofs);
    void print_pseudo(std::ofstream& ofs);
};

#endif // PSEUDO_H
