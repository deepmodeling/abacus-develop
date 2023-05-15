#include "module_basis/module_nao/numerical_radial.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include "module_base/constants.h"
#include "module_base/spherical_bessel_transformer.h"

using ModuleBase::PI;

NumericalRadial::NumericalRadial() {
    if (use_internal_transformer_) {
        sbt_ = new ModuleBase::SphericalBesselTransformer;
    }
}

NumericalRadial::NumericalRadial(NumericalRadial const& numerad) {
    this->symbol_ = numerad.symbol_;
    this->itype_ = numerad.itype_;
    this->ichi_ = numerad.ichi_;
    this->l_ = numerad.l_;

    this->nr_ = numerad.nr_;
    this->nk_ = numerad.nk_;

    this->is_fft_compliant_ = numerad.is_fft_compliant_;

    this->pr_ = numerad.pr_;
    this->pk_ = numerad.pk_;

    this->use_internal_transformer_ = numerad.use_internal_transformer_;

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

void NumericalRadial::build(
            const int l,                  
            const char r_or_k,            
            const int ngrid,              
            const double* const grid,     
            const double* const value,    
            const int p,
            const int itype,
            const int ichi,           
            const std::string symbol ) 
{
    assert(r_or_k == 'r' || r_or_k == 'k');
    assert( l >= 0 );
    assert( ngrid > 1 );
    assert( grid && value );

    symbol_ = symbol;
    itype_ = itype;
    ichi_= ichi;
    l_ = l;

    delete[] rgrid_;
    delete[] kgrid_;
    delete[] rvalue_;
    delete[] kvalue_;

    if (r_or_k == 'r') {
        nr_ = ngrid;
        pr_ = p;

        rgrid_ = new double[nr_];
        rvalue_ = new double[nr_];
        for (int ir = 0; ir != nr_; ++ir) {
            rgrid_[ir] = grid[ir];
            rvalue_[ir] = value[ir];
        }
    } else {
        nk_ = ngrid;
        pk_ = p;
        kgrid_ = new double[nk_];
        kvalue_ = new double[nk_];
        for (int ik = 0; ik != nk_; ++ik) {
            kgrid_[ik] = grid[ik];
            kvalue_[ik] = value[ik];
        }
    }
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
        // if no external transformer is provided
        if (!use_internal_transformer_) {
            sbt_ = new ModuleBase::SphericalBesselTransformer;
            use_internal_transformer_ = true;
        }
        // do nothing if an internal one is already in use
    }

    switch (update) {
        case 0: break;
        case 1: transform(true); break;
        case -1: transform(false); break;
        default: /* not supposed to happen */ ;
    }
}

void NumericalRadial::set_grid(
        const char r_or_k,       
        const int ngrid,         
        const double* const grid,
        const char mode   )
{
    assert( r_or_k == 'r' || r_or_k == 'k' );
    assert( mode == 'i' || mode == 't' );

    if (mode == 't') {
        if (r_or_k == 'r') {
            assert(kgrid_ && kvalue_);
            delete[] rgrid_;
            delete[] rvalue_;
            rgrid_ = new double[ngrid];
            rvalue_ = new double[ngrid];
            nr_ = ngrid;
            for (int ir = 0; ir != nr_; ++ir) {
                rgrid_[ir] = grid[ir];
            }
            check_fft_compliancy();
            transform(false);
        } else {
            assert(rgrid_ && rvalue_);
            delete[] kgrid_;
            delete[] kvalue_;
            kgrid_ = new double[ngrid];
            kvalue_ = new double[ngrid];
            nk_ = ngrid;
            for (int ik = 0; ik != nk_; ++ik) {
                kgrid_[ik] = grid[ik];
            }
            check_fft_compliancy();
            transform(true);
        }
    } else { // mode == 'i': interpolates the existing grid
        double* grid_tmp = new double[ngrid];
        double* value_tmp = new double[ngrid];
        for (int i = 0; i != ngrid; ++i) {
            grid_tmp[i] = grid[i];
        }

        if (r_or_k == 'r') {

            // TBD: interpolates rgrid_ & rvalues_ on grid to get value_tmp

            delete[] rgrid_;
            delete[] rvalue_;
            nr_ = ngrid;
            rgrid_ = std::move(grid_tmp);
            rvalue_ = std::move(value_tmp);
            check_fft_compliancy();
            transform(true);
        } else {

            // TBD: interpolates rgrid_ & rvalues_ on grid to get value_tmp

            delete[] kgrid_;
            delete[] kvalue_;
            nk_ = ngrid;
            kgrid_ = std::move(grid_tmp);
            kvalue_ = std::move(value_tmp);
            check_fft_compliancy();
            transform(false);
        }
    }
}

void NumericalRadial::set_value(const char r_or_k, const double* const value, const int p) {
    assert( r_or_k == 'r' || r_or_k == 'k' );
    if (r_or_k == 'r') {
        for (int ir = 0; ir != nr_; ++ir) {
            rvalue_[ir] = value[ir];
        }
        transform(true);
    } else {
        for (int ik = 0; ik != nk_; ++ik) {
            kvalue_[ik] = value[ik];
        }
        transform(false);
    }
}


void NumericalRadial::wipe(const char r_or_k) {

    assert(r_or_k == 'r' || r_or_k == 'k');

    if (r_or_k == 'r') {
        delete[] rgrid_;
        delete[] rvalue_;
        nr_ = 0;
        pr_ = 0;
    } else {
        delete[] kgrid_;
        delete[] rvalue_;
        nk_ = 0;
        pk_ = 0;
    }

    is_fft_compliant_ = false;
}


void NumericalRadial::save(std::string file) const {
    if (file.empty()) {
        file = symbol_ + "-" + std::to_string(l_) + "-" + std::to_string(ichi_) + ".dat";
    }
    std::ofstream ofs(file);
    if (ofs.is_open()) {
        ofs << "symbol " << symbol_ << std::endl;
        ofs << "l " << l_ << std::endl;
        ofs << "itype " << itype_ << std::endl;
        ofs << "ichi " << ichi_ << std::endl;
        ofs << "nr " << nr_ << std::endl;
        ofs << "nk " << nk_ << std::endl;
        ofs << "pr " << pr_ << std::endl;
        ofs << "pk " << pk_ << std::endl;
        ofs << "is_fft_compliant " << is_fft_compliant_ << std::endl;

        if (rgrid_ && rvalue_) {
            ofs << "rgrid " << "rvalue ";
            for (int ir = 0; ir != nr_; ++ir) {
                ofs << std::setw(20) << std::setprecision(15) 
                    << rgrid_[ir] << " " 
                    << rvalue_[ir] << std::endl;
            }
        }

        ofs << std::endl;

        if (kgrid_ && kvalue_) {
            ofs << "kgrid " << "kvalue ";
            for (int ik = 0; ik != nk_; ++ik) {
                ofs << std::setw(20) << std::setprecision(15) 
                    << kgrid_[ik] << " " << kvalue_[ik] << std::endl;
            }
        }
    }

    ofs.close();

}

void NumericalRadial::radtab(
            const char op,
            const NumericalRadial& ket,
            const int l, 
            double* const table,
            const bool deriv
) {
    assert(op == 'S' || op == 'I' || op == 'T' || op == 'U');
    assert(l >= 0);

    // currently only FFT-compliant grids are supported!
    // FFT-based transform requires that two NumericalRadial objects have exactly the same grid
    assert(this->is_fft_compliant_ && ket.is_fft_compliant_);
    assert(this->nr_ == ket.nr_);
    assert(this->rcut() == ket.rcut());

    double* ktmp = new double[nk_];
    for (int ik = 0; ik != nk_; ++ik) {
        ktmp[ik] = this->kvalue_[ik] * ket.kvalue_[ik];
    }

    int op_pk = 0;
    switch (op) {
        case 'T': op_pk = -2; break;
        case 'U': op_pk = 2; break;
    }

    if (deriv) { // derivative of radial table
        if (l == 0) {
            // j'_0(x) = -j_1(x)
            sbt_->radrfft(1, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk - 1);
            for (int ir = 0; ir != nr_; ++ir) {
                table[ir] *= -1;
            }
        } else {
            double* rtmp = new double[nr_];
            sbt_->radrfft(l+1, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk - 1);
            sbt_->radrfft(l-1, nk_, this->kcut(), ktmp, rtmp, this->pk_ + ket.pk_ + op_pk - 1);
            for (int ir = 0; ir != nr_; ++ir) {
                table[ir] = (l*rtmp[ir] - (l+1)*table[ir]) / (2*l+1);
            }
            delete[] rtmp;
        }
    } else {
        sbt_->radrfft(l, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk);
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
        pk_ = 0;
    } else {
        if (!rvalue_) {
            rvalue_ = new double[nr_];
        }
        sbt_->radrfft(l_, nk_, kgrid_[nk_-1], kvalue_, rvalue_, pk_);
        pr_ = 0;
    }
}

void NumericalRadial::check_fft_compliancy() {
    is_fft_compliant_ = false;

    if (!rgrid_ || !kgrid_ ) {
        return;
    }

    if (nr_ != nk_ ) {
        return;
    }

    if (nr_ < 2) {
        return;
    }

    double tol = 4.0*std::numeric_limits<double>::epsilon();

    double dr = rgrid_[nr_-1] / (nr_-1);
    for (int ir = 0; ir != nr_; ++ir) {
        if ( std::abs(ir*dr-rgrid_[ir]) > tol ) {
            return;
        }
    }

    double dk = kgrid_[nk_-1] / (nk_-1);
    if (std::abs(dr*dk-PI/(nr_-1)) > tol) {
        return;
    }

    for (int ik = 0; ik != nk_; ++ik) {
        if ( std::abs(ik*dk-kgrid_[ik]) > tol ) {
            return;
        }
    }

    is_fft_compliant_ = true;
}



