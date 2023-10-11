#ifndef PSEUDO_H
#define PSEUDO_H

#include "read_pp.h"
#include "../module_io/output.h"

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
      bool has_so;         // if .true. includes spin-orbit
      int nv;              // UPF file version number
      std::string psd;     // Element label
      std::string pp_type; // Pseudo type ( NC or US )
      bool tvanp;          // .true. if Ultrasoft
      bool nlcc;           // Non linear core corrections(bool)
      std::string xc_func; // Exch-Corr type
      int zv;              // z valence
      double etotps;       // total energy
      double ecutwfc;      // suggested cut-off for wfc
      double ecutrho;      // suggested cut-off for rho
      int lmax;            // maximum angular momentum component
      int mesh;            // number of point in the radial mesh
      int nchi;            // nwfc,number of wavefunctions
      int nbeta;           // number of projectors
      int nqlc;            // number of angular momenta in Q
      int kkbeta;          // kkbeta, point where the beta are zero

      std::string* els = nullptr; // els[nchi]
      int* lchi = nullptr;        // lchi[nchi]
      double* oc = nullptr;       // oc[nchi]

      double* jjj = nullptr;  // total angual momentum, jjj[nbeta]
      double* jchi = nullptr; // jchi(nwfc), added by zhengdy-soc
      int* nn = nullptr;

      // Local pseudopotentials
      double* vloc_at = nullptr; // [mesh], local potential( = pseudopot_upf.vloc )

      // <PP_MESH>
      double* r = nullptr;   // radial logaritmic mesh, r[0:mesh-1]
      double* rab = nullptr; // derivative of the radial mesh, rab[0:mesh-1]

      //<PP_NLCC>
      double* rho_atc = nullptr; // radial core charge density, rho_atc[0:mesh-1]

      //<PP_RHOATOM>
      double* rho_at = nullptr; // radial atomic charge density, rho_at[0:mesh-1]

      // <PP_PSWFC>
      ModuleBase::matrix chi; // radial atomic orbitals, chi(nchi, mesh)

      // other
      int msh;     // number of points up to rcut
      double rcut; // cut-off radius

      // <PP_BETA>
      int* lll = nullptr; // lll(nbeta), angular momentum of the beta function

      // <PP_DIJ>
      ModuleBase::matrix dion;  // dion(nbeta,nbeta)
      ModuleBase::matrix betar; // (nbeta, mesh), radial beta_{mu} functions

      // other
      int nh; // number of beta functions per atomic type

      // uspp
      ModuleBase::realArray qfuncl; // qfuncl(2*lmax+1,nbeta*(nbeta+1)/2,mesh) Q_{mu,nu}(|r|) function for |r|> r_L
      ModuleBase::matrix qqq;       // qqq(nbeta,nbeta) q_{mu,nu}

      void set_pseudo_h(const Pseudopot_upf& upf);
      void set_pseudo_atom(const Pseudopot_upf& upf);
      void set_pseudo_vl(const Pseudopot_upf& upf);
      void set_pseudo(const Pseudopot_upf& upf);

      void print_pseudo_h(std::ofstream& ofs);
      void print_pseudo_atom(std::ofstream& ofs);
      void print_pseudo_vl(std::ofstream& ofs);
      void print_pseudo(std::ofstream& ofs);
};

#endif // PSEUDO_H
