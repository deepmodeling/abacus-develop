#ifndef GRID_PARTITION_H
#define GRID_PARTITION_H

namespace Grid {
namespace Partition {

enum class Type {
    Becke,
    Stratmann,
    StratmannMod
};


/**
 * @brief Becke's partition weight.
 *
 * This function computes the normalized Becke's partition weight
 * for a grid point associated with a selected set of centers, given
 * the grid point's distance to centers and inter-center distance.
 *
 * @param   nR0     Total number of centers given by drR & dRR.
 * @param   drR     Distance between the grid point and centers.
 * @param   dRR     Distance between centers. dRR[I*nR0 + J] is the
 *                  distance between center I and J.
 * @param   nR      Number of centers involved in the weight calculation.
 *                  nR <= nR0. Length of iR.
 * @param   iR      Indices of centers involved.
 *                  Each element is a distinctive index in [0, nR0).
 * @param   c       iR[c] is the index of the center whom this grid point
 *                  belongs to.
 *
 * Reference:
 * Becke, A. D. (1988).
 * A multicenter numerical integration scheme for polyatomic molecules.
 * The Journal of chemical physics, 88(4), 2547-2553.
 *
 */
double w_becke(
    int nR0,
    double* drR,
    double* dRR,
    int nR,
    int* iR,
    int c
);

// Becke's cell function (iterated polynomial)
double s_becke(double mu);


/**
 * @brief Becke's partition weight with Stratmann's scheme.
 *
 * This function is similar to w_becke, but the cell function adopts
 * the one proposed by Stratmann et al, which enables some screening.
 *
 * @see w_becke
 *
 * @param   dRcRnn  Distance between the grid center and its nearest neighbor.
 *
 * Reference:
 * Stratmann, R. E., Scuseria, G. E., & Frisch, M. J. (1996).
 * Achieving linear scaling in exchange-correlation density functional
 * quadratures.
 * Chemical physics letters, 257(3-4), 213-223.
 *
 */
double w_stratmann(
    int nR0,
    double* drR,
    double* dRR,
    int nR,
    int* iR,
    int c,
    double dRcRnn,
    double a = 0.64
);

// Stratmann's piecewise cell function
double s_stratmann(double mu, double a);


/**
 * Knuth, F., Carbogno, C., Atalla, V., Blum, V., & Scheffler, M. (2015).
 * All-electron formalism for total energy strain derivatives and stress
 * tensor components for numeric atom-centered orbitals.
 * Computer Physics Communications, 190, 33-50.
 */
double s_stratmann_mod(double mu, double y, double a = 0.64, double b = 0.8);
double u_stratmann_mod(double y, double b); // used by s_stratmann_mod


} // end of namespace Partition
} // end of namespace Grid

#endif
