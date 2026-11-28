#include "source_estate/occ_mixer.h"

#include <algorithm>
#include <stdexcept>

void OccMatMixer::init(const UnitCell* cell,
                       const std::vector<int>* l_channel,
                       const std::vector<int>* flat_index,
                       const int nspin,
                       const int total_size)
{
    this->cell_ = cell;
    this->l_channel_ = l_channel;
    this->index_ = flat_index;
    this->nspin_ = nspin;
    this->uom_.resize(total_size, 0.0);
    this->uom_save_.resize(total_size, 0.0);
    this->raw_.resize(total_size, 0.0);
    this->initialized_.assign(total_size, 0);
    this->mixed_.assign(total_size, 0);
}

void OccMatMixer::seed_save(const OccupationMatrix& occmat)
{
    occmat.write_to_flat(*this->cell_, *this->l_channel_,
                          *this->index_, this->uom_save_);
    std::fill(this->initialized_.begin(), this->initialized_.end(), 1);
}

void OccMatMixer::begin_iter(OccupationMatrix& occmat)
{
    occmat.write_save_to_flat(*this->cell_, *this->l_channel_,
                              *this->index_, this->uom_save_);
}

void OccMatMixer::collect(const OccupationMatrix& occmat)
{
    occmat.write_to_flat(*this->cell_, *this->l_channel_,
                         *this->index_, this->uom_);
    this->raw_ = this->uom_;
    if (std::find(this->initialized_.begin(), this->initialized_.end(), 1) == this->initialized_.end())
    {
        this->uom_save_ = this->uom_;
    }
    std::fill(this->initialized_.begin(), this->initialized_.end(), 1);
}

void OccMatMixer::write_back(OccupationMatrix& occmat)
{
    occmat.read_from_flat(*this->cell_, *this->l_channel_,
                          *this->index_, this->uom_);
}

void OccMatMixer::mix_plain(OccupationMatrix& occmat, const double beta)
{
    elecstate::mix_occ_with_save(occmat.data(), occmat.data_save(),
                                 *this->cell_, *this->l_channel_,
                                 this->nspin_, beta);
}

int OccMatMixer::block_offset(const int atom, const int spin) const
{
    int offset = this->index_->at(atom);
    if (this->nspin_ == 2 && spin == 1)
    {
        offset += this->uom_.size() / 2;
    }
    return offset;
}

void OccMatMixer::mix_onsite_block(const int atom, const int spin,
                                  std::vector<double>& occupation, const double beta)
{
    const int offset = this->block_offset(atom, spin);
    if (offset < 0 || offset + occupation.size() > this->uom_.size())
    {
        throw std::out_of_range("Onsite occupation block exceeds mixing buffer");
    }
    const bool initialized = this->initialized_.at(offset) != 0;
    for (std::size_t i = 0; i < occupation.size(); ++i)
    {
        const int index = offset + i;
        this->raw_[index] = occupation[i];
        if (initialized)
        {
            occupation[i] = beta * occupation[i] + (1.0 - beta) * this->uom_save_[index];
        }
        this->uom_[index] = occupation[i];
        this->uom_save_[index] = occupation[i];
        this->initialized_[index] = 1;
        this->mixed_[index] = initialized;
    }
}

std::vector<double> OccMatMixer::raw_block(const int atom, const int spin, const int size) const
{
    const int offset = this->block_offset(atom, spin);
    if (this->mixed_.at(offset) == 0)
    {
        return {};
    }
    return std::vector<double>(this->raw_.begin() + offset, this->raw_.begin() + offset + size);
}

std::vector<double> OccMatMixer::mixed_block(const int atom, const int spin, const int size) const
{
    const int offset = this->block_offset(atom, spin);
    if (this->mixed_.at(offset) == 0)
    {
        return {};
    }
    return std::vector<double>(this->uom_.begin() + offset, this->uom_.begin() + offset + size);
}

void OccMatMixer::record_mixing(const bool applied)
{
    std::fill(this->mixed_.begin(), this->mixed_.end(), applied);
}
