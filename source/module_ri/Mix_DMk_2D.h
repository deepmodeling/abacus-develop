//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DMK_2D_H
#define MIX_DMK_2D_H

#include "Mix_Data.h"
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
	 * @param mixing_mode Mixing mode.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_mixing_mode(const Mixing_Mode mixing_mode);

	/**
	 * @brief Sets the mixing beta.
	 * @param mixing_beta Mixing beta.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_mixing_beta(const double mixing_beta);

	/**
	 * @brief Sets the pulay mixing coefficients from class Charge_Mixing.
	 * @tparam ChgMix Type of the charge mixing coefficient.
	 * @param iter Iteration step, start from 1.
	 * @param chr_mix Object of Charge_Mixing.
	 * @return Reference to the current object.
	 */
	template<typename ChgMix>
	Mix_DMk_2D &set_coef_pulay(const int iter, const ChgMix& chr_mix);

	/**
	 * @brief Mixes the double density matrix.
	 * @param dm Double Density matrix.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
	void mix(const std::vector<ModuleBase::matrix> &dm, const bool flag_restart);

	/**
	 * @brief Mixes the complex density matrix.
	 * @param dm Complex density matrix.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
	void mix(const std::vector<ModuleBase::ComplexMatrix> &dm, const bool flag_restart);

	/**
	 * @brief Returns the double density matrix.
	 * @return Double density matrices for each k-points.
	 */
	std::vector<const ModuleBase::matrix*> get_DMk_gamma_out() const;
	/**
	 * @brief Returns the complex density matrix.
	 * @return Complex density matrices for each k-points.
	 */
	std::vector<const ModuleBase::ComplexMatrix*> get_DMk_k_out() const;

private:
	std::vector<Mix_Data<ModuleBase::matrix>> mix_DMk_gamma;
	std::vector<Mix_Data<ModuleBase::ComplexMatrix>> mix_DMk_k;
	bool gamma_only;
};

#include "Mix_DMk_2D.hpp"

#endif