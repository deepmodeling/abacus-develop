#include "source_estate/occ_mixer.h"

void OccMatMixer::init(const UnitCell* cell,
                       const std::vector<int>* orbital_corr,
                       const std::vector<int>* flat_index,
                       const int nspin,
                       const int total_size)
{
    this->cell_ = cell;
    this->orbital_corr_ = orbital_corr;
    this->index_ = flat_index;
    this->nspin_ = nspin;
    this->uom_.resize(total_size, 0.0);
    this->uom_save_.resize(total_size, 0.0);
}

void OccMatMixer::seed_save(const OccupationMatrix& occmat)
{
    occmat.write_save_to_flat(*this->cell_, *this->orbital_corr_,
                              *this->index_, this->uom_save_);
}

void OccMatMixer::begin_iter(OccupationMatrix& occmat)
{
    // the caller has already snapshotted occ into occ_save via
    // OccupationMatrix::copy_to_save; here we only flatten that snapshot
    // into uom_save_ for the mixing history.
    occmat.write_save_to_flat(*this->cell_, *this->orbital_corr_,
                              *this->index_, this->uom_save_);
}

void OccMatMixer::collect(const OccupationMatrix& occmat)
{
    occmat.write_to_flat(*this->cell_, *this->orbital_corr_,
                         *this->index_, this->uom_);
}

void OccMatMixer::write_back(OccupationMatrix& occmat)
{
    occmat.read_from_flat(*this->cell_, *this->orbital_corr_,
                          *this->index_, this->uom_);
}

void OccMatMixer::mix_plain(OccupationMatrix& occmat, const double beta)
{
    elecstate::mix_occ_with_save(occmat.data(), occmat.data_save(),
                                 *this->cell_, *this->orbital_corr_,
                                 this->nspin_, beta);
}
