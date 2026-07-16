//=======================
// AUTHOR : Peize Lin
// DATE :   2023-05-09
//=======================

#include "Mix_DMk_2D.h"
#include "source_base/module_mixing/plain_mixing.h"
#include "source_base/tool_title.h"

#include <cassert>

template <typename Tdata>
Mix_DMk_2D<Tdata>::~Mix_DMk_2D<Tdata>()
{
    if(this->flag_del_mixing)
        delete this->mixing;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_nks(const int nks)
{
    this->mix_DMk.clear();
    this->mix_DMk.resize(nks);
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_mixing(Base_Mixing::Mixing* mixing_in)
{
    if(this->flag_del_mixing)
        delete this->mixing;
    this->mixing = mixing_in;
    this->flag_del_mixing = false;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::set_mixing_plain(const double& mixing_beta)
{
    if(this->flag_del_mixing)
        delete this->mixing;
    this->mixing = new Base_Mixing::Plain_Mixing(mixing_beta);
    this->flag_del_mixing = true;
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix(const std::vector<std::vector<Tdata>>& dm, const bool flag_restart)
{
    ModuleBase::TITLE("Mix_DMk_2D", "mix");
    if (flag_restart)
        { this->restart_all(dm); }
    else
        { this->mix_all(dm); }
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
void Mix_DMk_2D<Tdata>::restart_all(const std::vector<std::vector<Tdata>>& data_in)
{
    assert(this->mix_DMk.size() == data_in.size());
    assert(this->mixing != nullptr);
    for (int ik = 0; ik < data_in.size(); ++ik)
    {
        this->mix_DMk[ik].data_out = data_in[ik];
        this->mixing->init_mixing_data(this->mix_DMk[ik].mixing_data, data_in[ik].size(), sizeof(Tdata));
    }
}

template <typename Tdata>
void Mix_DMk_2D<Tdata>::mix_all(const std::vector<std::vector<Tdata>>& data_in)
{
    assert(this->mix_DMk.size() == data_in.size());
    assert(this->mixing != nullptr);
    for (int ik = 0; ik < data_in.size(); ++ik)
    {
        this->mixing->push_data(this->mix_DMk[ik].mixing_data, this->mix_DMk[ik].data_out.data(), data_in[ik].data(), nullptr, false);
        this->mixing->mix_data(this->mix_DMk[ik].mixing_data, this->mix_DMk[ik].data_out.data());
    }
}

template class Mix_DMk_2D<double>;
template class Mix_DMk_2D<std::complex<double>>;
