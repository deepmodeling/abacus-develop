#ifndef PSEUDO_NC_H
#define PSEUDO_NC_H

#include "read_pp.h"
#include "../module_io/output.h"

//-----------------------------------------
// read in norm conserving pseudopotentials
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
      std::string* els;    // els[nchi]
      int* lchi;           // lchi[nchi]
      double* oc;          // oc[nchi]

      double* jjj;  // total angual momentum, jjj[nbeta]
      double* jchi; // jchi(nwfc), added by zhengdy-soc
      int* nn;

      // Local pseudopotentials
      double* vloc_at; // [mesh], local potential( = pseudopot_upf.vloc )

      // <PP_MESH>
      double* r;   // radial logaritmic mesh, r[0:mesh-1]
      double* rab; // derivative of the radial mesh, rab[0:mesh-1]

      //<PP_NLCC>
      double* rho_atc; // radial core charge density, rho_atc[0:mesh-1]

      //<PP_RHOATOM>
      double* rho_at; // radial atomic charge density, rho_at[0:mesh-1]

      // <PP_PSWFC>
      ModuleBase::matrix chi; // radial atomic orbitals, chi(nchi, mesh)

      // other
      int msh;     // number of points up to rcut
      double rcut; // cut-off radius

      // <PP_BETA>
      int* lll;   // lll(nbeta), angular momentum of the beta function
      int kkbeta; // kkbeta(nbeta), point where the beta are zero

      // <PP_DIJ>
      ModuleBase::matrix dion;  // dion(nbeta,nbeta)
      ModuleBase::matrix betar; // (nbeta, mesh), radial beta_{mu} functions

      // other
      int nh; // number of beta functions per atomic type

      void set_pseudo_h(const Pseudopot_upf& upf);
      void set_pseudo_atom(const Pseudopot_upf& upf);
      void set_pseudo_vl(const Pseudopot_upf& upf);
      void set_pseudo_nc(const Pseudopot_upf& upf);

      void print_pseudo_h(std::ofstream& ofs);
      void print_pseudo_atom(std::ofstream& ofs);
      void print_pseudo_vl(std::ofstream& ofs);
      void print_pseudo_nc(std::ofstream& ofs);

};

#endif // PSEUDO_NC_H
