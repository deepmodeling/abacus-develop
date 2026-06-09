//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#include "Mix_DMk_2D.h"
#include "source_base/module_mixing/plain_mixing.h"
#include "source_base/tool_title.h"

#include <cassert>

template <typename Tdata>
Mix_DMk_2D<Tdata>& Mix_DMk_2D<Tdata>::set_nks(const int nks)
{
    ModuleBase::TITLE("Mix_DMk_2D", "set_nks");
    this->mix_DMk.clear();
    this->mix_DMk.resize(nks);
    return *this;
}

template <typename Tdata>
Mix_DMk_2D<Tdata>& Mix_DMk_2D<Tdata>::set_mixing(Base_Mixing::Mixing* mixing_in)
{
    ModuleBase::TITLE("Mix_DMk_2D", "set_mixing");
    this->mixing = mixing_in;
    this->separate_loop = (mixing_in == nullptr);
    return *this;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix(const std::vector<std::vector<Tdata>>& dm, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix");
    if (flag_restart)
    {
        this->restart_all(this->mix_DMk, dm);
    }
    else
    {
        this->mix_all(this->mix_DMk, dm);
    }
}

template <typename Tdata>
std::vector<const std::vector<Tdata>*> Mix_DMk_2D<Tdata>::get_DMk_out() const
{
    std::vector<const std::vector<Tdata>*> DMk_out(this->mix_DMk.size());
    for (int ik = 0; ik < this->mix_DMk.size(); ++ik)
        { DMk_out[ik] = &this->mix_DMk[ik].data_out; }
    return DMk_out;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::restart_one(typename Mix_DMk_2D<Tdata>::DMk_Mix_Data& data,
                                    const typename Mix_DMk_2D<Tdata>::Tmatrix& data_in,
                                    Base_Mixing::Mixing& mixing)
{
    data.data_out = data_in;
    const int length = static_cast<int>(data_in.size());
    mixing.init_mixing_data(data.mixing_data, length, sizeof(Tdata));
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix_one(typename Mix_DMk_2D<Tdata>::DMk_Mix_Data& data,
                                const typename Mix_DMk_2D<Tdata>::Tmatrix& data_in,
                                Base_Mixing::Mixing& mixing)
{
    mixing.push_data(data.mixing_data, data.data_out.data(), data_in.data(), nullptr, false);
    mixing.mix_data(data.mixing_data, data.data_out.data());
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::restart_all(std::vector<typename Mix_DMk_2D<Tdata>::DMk_Mix_Data>& data_out,
                                    const std::vector<typename Mix_DMk_2D<Tdata>::Tmatrix>& data_in)
{
    assert(data_out.size() == data_in.size());
    if (this->separate_loop)
    {
        Base_Mixing::Plain_Mixing plain_mixing(1.0);
        for (int ik = 0; ik < data_in.size(); ++ik)
        {
            restart_one(data_out[ik], data_in[ik], plain_mixing);
        }
    }
    else
    {
        assert(this->mixing != nullptr);
        for (int ik = 0; ik < data_in.size(); ++ik)
        {
            restart_one(data_out[ik], data_in[ik], *this->mixing);
        }
    }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix_all(std::vector<typename Mix_DMk_2D<Tdata>::DMk_Mix_Data>& data_out,
                                const std::vector<typename Mix_DMk_2D<Tdata>::Tmatrix>& data_in)
{
    assert(data_out.size() == data_in.size());
    if (this->separate_loop)
    {
        Base_Mixing::Plain_Mixing plain_mixing(1.0);
        for (int ik = 0; ik < data_in.size(); ++ik)
        {
            mix_one(data_out[ik], data_in[ik], plain_mixing);
        }
    }
    else
    {
        assert(this->mixing != nullptr);
        for (int ik = 0; ik < data_in.size(); ++ik)
        {
            mix_one(data_out[ik], data_in[ik], *this->mixing);
        }
    }
}

template class Mix_DMk_2D<double>;
template class Mix_DMk_2D<std::complex<double>>;
