#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_TABLE_PHI_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_BASIS_MODULE_AO_ORB_TABLE_PHI_H

#include "ORB_atomic_lm.h"
#include "ORB_read.h"
#include "module_base/intarray.h"
#include "module_base/sph_bessel_recursive.h"

#include <set>

class ORB_table_phi
{
  public:
    ORB_table_phi();
    ~ORB_table_phi();

    void allocate(const int& ntype,      ///< number of atom types
                  const int& lmax_in,    ///< max L used to calculate overlap
                  const int& kmesh_in,   ///< kpoints, for integration in k space
                  const double& Rmax_in, ///< max value of radial table
                  const double& dR_in,   ///< delta R, for making radial table
                  const double& dk_in    ///< delta k, for integration in k space
    );

    void init_Table(LCAO_Orbitals& orb);

    void Destroy_Table(LCAO_Orbitals& orb);

    /**
     * Five dimension:
     *-----------------------------
     * (1) 0: normal (S(R)) ; 1: derivative( dS/dR )
     *
     * (2) pairs type number.
     *
     * (3) pairs chi
     *
     * (4) Max angular momentum: L.
     *
     * (5) Distance between atoms: R.
     */
    double***** Table_SR;
    double***** Table_TR;

    bool overlap_table_allocated;
    bool kinetic_table_allocated;

    // Wenfei 2021-8-26, plot table elements against R
    void plot_table(const std::string filename, const int rmesh, double* column);

    // Peize Lin add 2017-04-24, and change all jlx in this class
    ModuleBase::Sph_Bessel_Recursive::D2* pSB = nullptr;

    ///
    /// make the index, in order to get the element from Table_SR and Table_TR quickly.
    ///

    ///
    /// OV stands for 'overlap'
    ///
    /// T stands for atom type.
    ///
    void init_OV_Tpair(LCAO_Orbitals& orb);
    ///
    /// O stands for orbitals.
    ///
    void init_OV_Opair(LCAO_Orbitals& orb);

    int OV_nTpairs;
    ModuleBase::IntArray OV_Tpair;
    ModuleBase::IntArray OV_Opair;
    ModuleBase::IntArray OV_L2plus1;

    double dr;
    int Rmesh;

  private:
    // variables
    int ntype;

    int lmax;

    double Rmax;

    double dk;

    int nlm;

    int kmesh;

    double* kpoint;

    double* r;

    double* rab;

    double* kab;

    //------------------------------
    // sizes of table
    int nelem_;                 // number of elements
    std::vector<int> lmax_;     // lmax of each element
    std::vector<int> nchi_tot_; // total nchi of each element

    // automatically deallocate Table_DSR using lmax_d_, lmax_ & nchi_pairs_
    // called by destructor
    void _destroy_table();
};
#endif
