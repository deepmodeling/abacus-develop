#ifndef PSEUDOPOT_UPF_H
#define PSEUDOPOT_UPF_H

#include <string>

#include "atom_pseudo.h"
#include "module_base/matrix.h"
#include "module_base/realarray.h"

class Pseudopot_upf
{
public:
	//PP_INFO
	//PP_HEADER
	//PP_MESH
	//PP_NLCC
	//PP_LOCAL
	//PP_NONLOCAL
	//PP_PSWFC
	//PP_PSRHOATOM
	//addinfo

	Pseudopot_upf();
	~Pseudopot_upf();

    std::string relativistic; // relativistic: no, scalar, full
    int lmax_rho;             // maximum angular momentum component in rho (should be 2*lmax)
    double xmin;              // the minimum x of the linear mesh
    double rmax;              // the maximum radius of the mesh
    double zmesh;             // the nuclear charge used for mesh
    double dx;                // the deltax of the linear mesh
                              // The radial grid is: r(i+1) = exp(xmin+i*dx)/zed  a.u.
    int lloc;                 // L of channel used to generate local potential
                              //   (if < 0 it was generated by smoothing AE potential)
    // double rcloc;             // vloc = v_ae for r > rcloc
    bool q_with_l;            // if .true. qfunc is pseudized in
    int nqf;                  // number of Q coefficients
    // bool has_wfc;             // if true, UPF contain AE and PS wfc for each beta

    // need 'new' and 'delete'
    bool coulomb_potential = false;  // coulomb potentail : z/r
    ModuleBase::matrix chi;          // chi(nwfc,mesh) atomic wavefcts
    std::vector<int> kbeta = {};            // kbeta(nbeta):number of mesh points for projector i (must be .le. mesh )
    std::vector<std::string> els_beta = {}; // els_beta(nwfc):label for the beta
    std::vector<int> nchi = {};             // nchi(nwfc) value of pseudo-n for wavefcts
    std::vector<double> epseu = {};         // epseu(nwfc) pseudo one-particle energy
    std::vector<double> rcut_chi = {};      // rcut_chi(nwfc) cutoff inner radius
    std::vector<double> rcutus_chi = {};    // rcutus_chi(nwfc) ultrasoft outer radius
    std::vector<double> rinner = {};        // rinner(2*lmax+1) r_L
    ModuleBase::matrix qfunc;        // qfunc(nbeta*(nbeta+1)/2,mesh) Q_{mu,nu}(|r|) function for |r|> r_L
    ModuleBase::realArray qfcoef;    // qfcoef(nbeta,nbeta,2*lmax+1,nqf) coefficients for Q for |r|<r_L
    // ModuleBase::matrix aewfc;        // wfc(nbeta,mesh) all-electron wfc
    // ModuleBase::matrix pswfc;        // wfc(nbeta,mesh) pseudo wfc
    std::vector<double> rcut = {};          // cut-off radius(nbeta)
    std::vector<double> rcutus = {};        // ultrasoft cut-off radius (nbeta)

    int nd; // nl_5 // Number of nonzero Dij

    // the followings are for the vwr format
    int spd_loc;
    int iTB_s;
    int iTB_p;
    int iTB_d;

    // return error
    int init_pseudo_reader(const std::string& fn, std::string& type, Atom_pseudo& pp);
    void print_pseudo_upf(std::ofstream& ofs, Atom_pseudo& pp);

    int average_p(const double& lambda, Atom_pseudo& pp); // zhengdy add 2020-10-20
    void set_empty_element(Atom_pseudo& pp);            // Peize Lin add for bsse 2022.04.07
    void set_upf_q(Atom_pseudo& pp);                    // liuyu add 2023-09-21
    void complete_default(Atom_pseudo& pp);

  private:
    int set_pseudo_type(const std::string& fn, std::string& type);
    std::string& trim(std::string& in_str);
    std::string trimend(std::string& in_str);

    int read_pseudo_upf(std::ifstream& ifs, Atom_pseudo& pp);
    int read_pseudo_vwr(std::ifstream& ifs, Atom_pseudo& pp);
    int read_pseudo_blps(std::ifstream& ifs, Atom_pseudo& pp); // sunliang added 2021.07.08
    void read_pseudo_header(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_mesh(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_nlcc(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_local(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_nl(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_pswfc(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_rhoatom(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_addinfo(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_so(std::ifstream& ifs, Atom_pseudo& pp);

    // upf201
    int read_pseudo_upf201(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_upf201_header(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_upf201_mesh(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_upf201_nonlocal(std::ifstream& ifs, Atom_pseudo& pp);
    void read_pseudo_upf201_pswfc(std::ifstream& ifs, Atom_pseudo& pp);
    // void read_pseudo_upf201_fullwfc(std::ifstream& ifs);
    void read_pseudo_upf201_so(std::ifstream& ifs, Atom_pseudo& pp);
    void getnameval(std::ifstream&, int&, std::string*, std::string*);

    /**
     * @brief Computes the Q function from its polynomial expansion (r < rinner)
     * @param nqf number of polynomial coefficients
     * @param mesh number of mesh points
     * @param l angular momentum
     * @param n additional exponent, result is multiplied by r^n
     * @param qfcoef polynomial coefficients
     * @param r radial mesh
     * @param rho output: r^n * Q(r)
     */
    void setqfnew(const int& nqf,
                  const int& mesh,
                  const int& l,
                  const int& n,
                  const double* qfcoef,
                  const double* r,
                  double* rho);

    // complete default
    // void complete_default(Atom_pseudo& pp);
    void complete_default_h(Atom_pseudo& pp);
    void complete_default_atom(Atom_pseudo& pp);
    void complete_default_vl(Atom_pseudo& pp);
};

#endif //pseudopot_upf class
