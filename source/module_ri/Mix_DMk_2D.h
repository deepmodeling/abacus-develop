//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DMK_2D_H
#define MIX_DMK_2D_H

#include "Mix_Matrix.h"
#include "module_base/module_mixing/mixing.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

#include <vector>

class Mix_DMk_2D
{
public:
	/**
	 * @brief Sets the number of k-points and gamma_only flag.
	 * @param nks Number of k-points.
	 * @param gamma_only_in Flag indicating if only gamma point is mixed.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_nks(const int nks, const bool gamma_only_in);

	/**
	 * @brief Sets the mixing mode.
	 * @param Mixing Mixing pointer.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_mixing(Base_Mixing::Mixing* mixing_in);

	/**
	 * @brief Sets the mixing beta.
	 * @param mixing_beta Mixing beta.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_mixing_beta(const double mixing_beta);

	/**
	 * @brief Mixes the double density matrix.
	 * @param dm Double Density matrix.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
    void mix(const std::vector<std::vector<double>>& dm, const bool flag_restart);

	/**
	 * @brief Mixes the complex density matrix.
	 * @param dm Complex density matrix.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
    void mix(const std::vector<std::vector<std::complex<double>>>& dm, const bool flag_restart);

	/**
	 * @brief Returns the double density matrix.
	 * @return Double density matrices for each k-points.
	 */
    std::vector<const std::vector<double>*> get_DMk_gamma_out() const;
	/**
	 * @brief Returns the complex density matrix.
	 * @return Complex density matrices for each k-points.
	 */
    std::vector<const std::vector<std::complex<double>>*> get_DMk_k_out() const;

private:
    std::vector<Mix_Matrix<std::vector<double>>> mix_DMk_gamma;
    std::vector<Mix_Matrix<std::vector<std::complex<double>>>> mix_DMk_k;
	bool gamma_only;
};

#endif