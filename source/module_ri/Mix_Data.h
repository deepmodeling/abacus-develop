//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DATA_H
#define MIX_DATA_H

#include <deque>
#include <vector>

enum class Mixing_Mode{No, Plain, Pulay};

	class Charge_Mixing;

template<typename Tdata>
class Mix_Data
{
public:
	Mixing_Mode mixing_mode = Mixing_Mode::No;
	double mixing_beta = 1.0;
	const Tdata &get_data_out() const { return this->data_out; }

	/**
	 * @brief Mixes the input data according to the set mixing mode.
	 * @param data_in Input data to be mixed.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
	void mix(const Tdata &data_in, const bool flag_restart);

	/**
	 * @brief Sets the pulay mixing coefficients from class Charge_Mixing.
	 * @tparam ChgMix Type of the charge mixing coefficient.
	 * @param iter Iteration step, start from 1.
	 * @param chr_mix Object of Charge_Mixing.
	 */
    template<typename ChgMix>
    void set_coef_pulay(const int iter, const ChgMix& chr_mix);

private:
	Tdata data_out;
	std::deque<Tdata> data_pulay_list;
	std::vector<double> coef_pulay_list;

	void pulay_mixing(const Tdata &data_in, const bool flag_restart);
};

#include "Mix_Data.hpp"

#endif