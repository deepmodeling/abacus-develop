#include "module_basis//module_nao/numerical_radial.h"
#include "module_base/spherical_bessel_transformer.h"

NumericalRadial::NumericalRadial() {
    // internal SphericalBesselTransformer
    sbt_ = new ModuleBase::SphericalBesselTransformer;
}

NumericalRadial::NumericalRadial(NumericalRadial const& numerad) {
    this->l_ = numerad.l_;
    this->nr_ = numerad.nr_;
    this->nk_ = numerad.nk_;

    this->pr_ = numerad.pr_;
    this->pk_ = numerad.pk_;

    this->has_uniform_rgrid_ = numerad.has_uniform_rgrid_;
    this->has_uniform_kgrid_ = numerad.has_uniform_kgrid_;
    this->is_fft_compliant_ = numerad.is_fft_compliant_;

    // deep copy
    this->rgrid_ = new double[nr_];
    this->rvalue_= new double[nr_];
    for (int ir = 0; ir != nr_; ++ir) {
        this->rgrid_[ir] = numerad.rgrid_[ir];
        this->rvalue_[ir] = numerad.rvalue_[ir];
    }

    this->kgrid_ = new double[nk_];
    this->kvalue_= new double[nk_];
    for (int ik = 0; ik != nk_; ++ik) {
        this->kgrid_[ik] = numerad.kgrid_[ik];
        this->kvalue_[ik] = numerad.kvalue_[ik];
    }

    this->use_internal_transformer_ = numerad.use_internal_transformer_;
    if (use_internal_transformer_) {
        this->sbt_ = new ModuleBase::SphericalBesselTransformer;
    } else {
        this->sbt_ = numerad.sbt_;
    }

}

NumericalRadial::~NumericalRadial() {
    delete[] rgrid_;
    delete[] kgrid_;
    delete[] rvalue_;
    delete[] kvalue_;

    if (use_internal_transformer_) {
        delete sbt_;
    }
}

void NumericalRadial::wipe(const char r_or_k) {

    assert(r_or_k == 'r' || r_or_k == 'k');

    if (r_or_k == 'r') {
        delete[] rgrid_;
        delete[] rvalue_;
        nr_ = 0;
        pr_ = 0;
        has_uniform_rgrid_ = false;
    } else {
        delete[] kgrid_;
        delete[] rvalue_;
        nk_ = 0;
        pk_ = 0;
        has_uniform_kgrid_ = false;
    }

    is_fft_compliant_ = false;
}


void NumericalRadial::set_transformer(ModuleBase::SphericalBesselTransformer* sbt, int update) {
    assert( update == 0 || update == 1 || update == -1 );
    if (sbt) {
        //! if an external transformer is provided
        if (use_internal_transformer_) {
            delete sbt_;
            use_internal_transformer_ = false;
        }
        sbt_ = sbt;
    } else {
        // if no external transformer is provided (implies the use of an internal one)
        if (!use_internal_transformer_) {
            sbt_ = new ModuleBase::SphericalBesselTransformer;
            use_internal_transformer_ = true;
        }
        // do nothing if an internal one is used already
    }

    switch (update) {
        case 0: break;
        case 1: transform(true); break;
        case -1: transform(false); break;
        default: /* not supposed to happen */ ;
    }
}

void NumericalRadial::transform(const bool forward) {
    if (forward) {
        assert( rgrid_ && rvalue_ );
        if (!kgrid_) {
            return;
        }
    } else {
        assert( kgrid_ && kvalue_ );
        if (!rgrid_) {
            return;
        }
    }

    // currently we only support FFT-compliant grid!
    assert( is_fft_compliant_ );

    if (forward) {
        if (!kvalue_) {
            kvalue_ = new double[nk_];
        }
        sbt_->radrfft(l_, nr_, rgrid_[nr_-1], rvalue_, kvalue_, pr_);
    } else {
        sbt_->radrfft(l_, nk_, kgrid_[nk_-1], kvalue_, rvalue_, pk_);
        if (!rvalue_) {
            rvalue_ = new double[nr_];
        }
    }
}

