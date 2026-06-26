//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_DMK_2D_H
#define MIX_DMK_2D_H

#include "source_base/module_mixing/mixing.h"

#include <complex>
#include <vector>

template <typename Tdata>
class Mix_DMk_2D
{
public:
	~Mix_DMk_2D<Tdata>();

	/**
	 * @brief Sets the number of k-points.
	 * @param nks Number of k-points.
	 */
	void set_nks(const int nks);

	/**
	 * @brief Sets the mixing mode.
	 * @param Mixing Mixing pointer.
	 */
	void set_mixing(Base_Mixing::Mixing* mixing_in);

	/**
	 * @brief Sets Base_Mixing::Plain_Mixing.
	 * @param mixing_beta mixing beta for plain mixing.
	 */
	void set_mixing_plain(const double& mixing_beta);

	/**
	 * @brief Mixes the density matrix.
	 * @param dm Density matrix.
	 * @param flag_restart Flag indicating whether restart mixing.
	 */
    void mix(const std::vector<std::vector<Tdata>>& dm, const bool flag_restart);

	/**
	 * @brief Returns the density matrix.
	 * @return Density matrices for each k-points.
	 */
    std::vector<const std::vector<Tdata>*> get_DMk_out() const;

private:
    struct DMk_Mix_Data
    {
        std::vector<Tdata> data_out;
        Base_Mixing::Mixing_Data mixing_data;
    };

    void restart_all(const std::vector<std::vector<Tdata>>& data_in);

    void mix_all(const std::vector<std::vector<Tdata>>& data_in);

    std::vector<DMk_Mix_Data> mix_DMk;
    Base_Mixing::Mixing* mixing = nullptr;
	bool flag_del_mixing = false;
};

#endif
