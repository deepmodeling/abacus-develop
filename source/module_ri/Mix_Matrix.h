//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#ifndef MIX_MATRIX_H
#define MIX_MATRIX_H

#include <vector>

#include "module_base/module_mixing/mixing.h"
#include "module_base/module_mixing/plain_mixing.h"
template <typename Tdata>
class Mix_Matrix
{
  public:
    /**
     * @brief init Mix_Matrix with Mixing pointer
     *
     * @param mixing_in Mixing pointer
     */
    void init(Base_Mixing::Mixing* mixing_in)
    {
        this->mixing = mixing_in;
        if (this->mixing == nullptr)
            this->separate_loop = true;
    }
    /**
     * @brief Get the data out object
     *
     * @return const Tdata&
     */
    const Tdata& get_data_out() const
    {
        return this->data_out;
    }

    /**
     * @brief Mixes the input data according to the set mixing mode.
     * @param data_in Input data to be mixed.
     * @param flag_restart Flag indicating whether restart mixing.
     */
    void mix(const Tdata& data_in, const bool flag_restart);

    double mixing_beta = 1.0;

  private:
    Tdata data_out;
    Base_Mixing::Mixing* mixing = nullptr;
    Base_Mixing::Mixing_Data matrix_data;
    bool separate_loop = false;
};
#endif