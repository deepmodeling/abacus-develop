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
	/**
	 * @brief Sets the number of k-points.
	 * @param nks Number of k-points.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_nks(const int nks);

	/**
	 * @brief Sets the mixing mode.
	 * @param Mixing Mixing pointer.
	 * @return Reference to the current object.
	 */
	Mix_DMk_2D &set_mixing(Base_Mixing::Mixing* mixing_in);

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
    using Tmatrix = std::vector<Tdata>;

    struct DMk_Mix_Data
    {
        Tmatrix data_out;
        Base_Mixing::Mixing_Data mixing_data;
    };

    static void restart_one(DMk_Mix_Data& data,
                            const Tmatrix& data_in,
                            Base_Mixing::Mixing& mixing);

    static void mix_one(DMk_Mix_Data& data,
                        const Tmatrix& data_in,
                        Base_Mixing::Mixing& mixing);

    void restart_all(std::vector<DMk_Mix_Data>& data_out,
                     const std::vector<Tmatrix>& data_in);

    void mix_all(std::vector<DMk_Mix_Data>& data_out,
                 const std::vector<Tmatrix>& data_in);

    std::vector<DMk_Mix_Data> mix_DMk;
    Base_Mixing::Mixing* mixing = nullptr;
    bool separate_loop = false;
};

#endif
