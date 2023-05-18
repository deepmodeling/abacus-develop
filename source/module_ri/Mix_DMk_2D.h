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
    Mix_DMk_2D& set_mixing_beta(const double mixing_beta);
    template<typename ChgMix>
    Mix_DMk_2D& set_coef_pulay(const int iter, const ChgMix& chr_mix);

/**
     * @brief Mixes the density matrix.
     * @param dm double Density matrix.
     * @param flag_restart Flag indicating if it is a restart.
     */
	void mix(const std::vector<ModuleBase::matrix> &dm, const bool flag_restart);

/**
     * @brief Mixes the complex density matrix.
     * @param dm Complex density matrix.
     * @param flag_restart Flag indicating if it is a restart.
     */
	void mix(const std::vector<ModuleBase::ComplexMatrix> &dm, const bool flag_restart);

/**
     * @brief Returns the gamma density matrix.
     * @return Vector of pointers to gamma density matrices.
     */
	std::vector<const ModuleBase::matrix*> get_DMk_gamma_out() const;
/**
     * @brief Returns the k-point density matrix.
     * @return Vector of pointers to k-point density matrices.
     */
	std::vector<const ModuleBase::ComplexMatrix*> get_DMk_k_out() const;

private:
	std::vector<Mix_Data<ModuleBase::matrix>> mix_DMk_gamma;
    std::vector<Mix_Data<ModuleBase::ComplexMatrix>> mix_DMk_k;
    bool gamma_only;
};

template<typename ChgMix>
Mix_DMk_2D& Mix_DMk_2D::set_coef_pulay(const int iter, const ChgMix& chr_mix)
{
	ModuleBase::TITLE("Mix_DMk_2D","set_coef_pulay");
	if(gamma_only)
		for(Mix_Data<ModuleBase::matrix> &mix_one : this->mix_DMk_gamma)
			mix_one.set_coef_pulay(iter, chr_mix);
	else
		for(Mix_Data<ModuleBase::ComplexMatrix> &mix_one : this->mix_DMk_k)
			mix_one.set_coef_pulay(iter, chr_mix);
	return *this;
}

#endif