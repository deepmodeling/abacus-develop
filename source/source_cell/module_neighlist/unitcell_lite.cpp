#include "unitcell_lite.h"

#include <cassert>

// === AtomProvider interface implementation ===

double UnitCellLite::get_lat0() const {
    return lat0_;
}

double UnitCellLite::get_omega() const {
    return omega_;
}

const ModuleBase::Matrix3& UnitCellLite::get_latvec() const {
    return latvec_;
}

int UnitCellLite::get_natom() const {
    return nat_;
}

int UnitCellLite::get_na(int i) const {
    assert(i >= 0 && i < ntype_);
    return na_[i];
}

int UnitCellLite::get_ntype() const {
    return ntype_;
}

ModuleBase::Vector3<double> UnitCellLite::get_tau(int i, int j) const {
    assert(i >= 0 && i < ntype_);
    assert(j >= 0 && j < na_[i]);
    if (i == 0) {
        return tau_[j];
    }
    return tau_[naa_[i - 1] + j];
}

// === Setter methods ===

void UnitCellLite::set_lat0(double lat0) {
    lat0_ = lat0;
}

void UnitCellLite::set_omega(double omega) {
    omega_ = omega;
}

void UnitCellLite::set_latvec(const ModuleBase::Matrix3& latvec) {
    latvec_ = latvec;
}

void UnitCellLite::set_lattice(double lat0, double omega, const ModuleBase::Matrix3& latvec) {
    lat0_ = lat0;
    omega_ = omega;
    latvec_ = latvec;
}

void UnitCellLite::set_atoms(int ntype,
                              const std::vector<int>& na,
                              const std::vector<ModuleBase::Vector3<double>>& tau) {
    assert(ntype >= 0);
    assert(na.size() == static_cast<size_t>(ntype));

    ntype_ = ntype;
    na_ = na;
    tau_ = tau;

    // compute total number of atoms
    nat_ = 0;
    for (int i = 0; i < ntype_; ++i) {
        nat_ += na_[i];
    }
    assert(tau_.size() == static_cast<size_t>(nat_));

    // compute cumulative counts
    compute_naa_();
}

// === Internal methods ===

void UnitCellLite::compute_naa_() {
    naa_.resize(na_.size());
    if (naa_.size() > 0) {
        naa_[0] = na_[0];
    }
    for (size_t i = 1; i < naa_.size(); ++i) {
        naa_[i] = naa_[i - 1] + na_[i];
    }
}