#ifndef SEP_H
#define SEP_H

#include <fstream>
#include <string>

/**
 * Sep Potential for DFT-1/2 etc.
 *
 * Sep Potential
 */
class SepPot
{
  public:
    SepPot();
    ~SepPot();

    bool is_enable = false;
    double r_in = 0.0;      /**< cut-off radius inner */
    double r_out = 0.0;     /**< cut-off radius outter */
    double r_power = 20.0;  /**< shell function exp factor */
    double enhence_a = 1.0; /**< scale sep potential */
    std::string label;      /**< element nameof sep  */
    std::string xc_type;    /**< Exch-Corr type */
    std::string orbital;    /** atomic angular moment s,p,d,f */
    int mesh = 0;           /**< number of points in radial mesh */
    int strip_elec = 0;     /**< strip electron amount 1->0.01 50->0.5 */
    double* r = nullptr;    /**< ridial mesh */
    double* rv = nullptr;   /**< sep potential, but rV, unit: Ry */

    int read_sep(std::ifstream& is);
    void print_sep_info(std::ofstream& ofs);
    void print_sep_vsep(std::ofstream& ofs);
#ifdef __MPI
    void bcast_sep();
#endif /* ifdef __MPI */
};

#endif /* ifndef SEP_H */
