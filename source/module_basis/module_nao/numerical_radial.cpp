#include "module_basis/module_nao/numerical_radial.h"

#include <fstream>
#include <iomanip>
#include <limits>

#include "module_base/constants.h"
#include "module_base/mathzone_add1.h"
#include "module_base/spherical_bessel_transformer.h"

using ModuleBase::Mathzone_Add1;
using ModuleBase::PI;

NumericalRadial::NumericalRadial()
{
    if (use_internal_transformer_)
    {
        sbt_ = new ModuleBase::SphericalBesselTransformer;
    }
}

NumericalRadial::NumericalRadial(const NumericalRadial& other)
{
    this->symbol_ = other.symbol_;
    this->itype_ = other.itype_;
    this->ichi_ = other.ichi_;
    this->l_ = other.l_;

    this->nr_ = other.nr_;
    this->nk_ = other.nk_;

    this->is_fft_compliant_ = other.is_fft_compliant_;

    this->pr_ = other.pr_;
    this->pk_ = other.pk_;

    this->use_internal_transformer_ = other.use_internal_transformer_;

    // deep copy
    if (other.ptr_rgrid())
    {
        this->rgrid_ = new double[nr_];
        this->rvalue_ = new double[nr_];
        for (int ir = 0; ir != nr_; ++ir)
        {
            this->rgrid_[ir] = other.rgrid_[ir];
            this->rvalue_[ir] = other.rvalue_[ir];
        }
    }

    if (other.ptr_kgrid())
    {
        this->kgrid_ = new double[nk_];
        this->kvalue_ = new double[nk_];
        for (int ik = 0; ik != nk_; ++ik)
        {
            this->kgrid_[ik] = other.kgrid_[ik];
            this->kvalue_[ik] = other.kvalue_[ik];
        }
    }

    if (use_internal_transformer_)
    {
        this->sbt_ = new ModuleBase::SphericalBesselTransformer;
    }
    else
    {
        this->sbt_ = other.sbt_;
    }
}

NumericalRadial& NumericalRadial::operator=(const NumericalRadial& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    this->symbol_ = rhs.symbol_;
    this->itype_ = rhs.itype_;
    this->ichi_ = rhs.ichi_;
    this->l_ = rhs.l_;

    this->nr_ = rhs.nr_;
    this->nk_ = rhs.nk_;

    this->is_fft_compliant_ = rhs.is_fft_compliant_;

    this->pr_ = rhs.pr_;
    this->pk_ = rhs.pk_;

    this->use_internal_transformer_ = rhs.use_internal_transformer_;

    // deep copy
    if (rhs.ptr_rgrid())
    {
        this->rgrid_ = new double[nr_];
        this->rvalue_ = new double[nr_];
        for (int ir = 0; ir != nr_; ++ir)
        {
            this->rgrid_[ir] = rhs.rgrid_[ir];
            this->rvalue_[ir] = rhs.rvalue_[ir];
        }
    }

    if (rhs.ptr_kgrid())
    {
        this->kgrid_ = new double[nk_];
        this->kvalue_ = new double[nk_];
        for (int ik = 0; ik != nk_; ++ik)
        {
            this->kgrid_[ik] = rhs.kgrid_[ik];
            this->kvalue_[ik] = rhs.kvalue_[ik];
        }
    }

    if (use_internal_transformer_)
    {
        this->sbt_ = new ModuleBase::SphericalBesselTransformer;
    }
    else
    {
        this->sbt_ = rhs.sbt_;
    }

    return *this;
}

NumericalRadial::~NumericalRadial()
{
    delete[] rgrid_;
    delete[] kgrid_;
    delete[] rvalue_;
    delete[] kvalue_;

    if (use_internal_transformer_)
    {
        delete sbt_;
    }
}

void NumericalRadial::build(const int l,
                            const bool for_r_space,
                            const int ngrid,
                            const double* const grid,
                            const double* const value,
                            const int p,
                            const int itype,
                            const int ichi,
                            const std::string symbol)
{
    assert(l >= 0);
    assert(ngrid > 1);
    assert(grid && value);

    symbol_ = symbol;
    itype_ = itype;
    ichi_ = ichi;
    l_ = l;

    delete[] rgrid_;
    delete[] kgrid_;
    delete[] rvalue_;
    delete[] kvalue_;
    rgrid_ = nullptr;
    kgrid_ = nullptr;
    rvalue_ = nullptr;
    kvalue_ = nullptr;

    if (for_r_space)
    {
        nr_ = ngrid;
        pr_ = p;
        rgrid_ = new double[nr_];
        rvalue_ = new double[nr_];
        for (int ir = 0; ir != nr_; ++ir)
        {
            rgrid_[ir] = grid[ir];
            rvalue_[ir] = value[ir];
        }
    }
    else
    {
        nk_ = ngrid;
        pk_ = p;
        kgrid_ = new double[nk_];
        kvalue_ = new double[nk_];
        for (int ik = 0; ik != nk_; ++ik)
        {
            kgrid_[ik] = grid[ik];
            kvalue_[ik] = value[ik];
        }
    }
}

void NumericalRadial::set_transformer(ModuleBase::SphericalBesselTransformer* sbt, int update)
{

    assert(update == 0 || update == 1 || update == -1);

    if (sbt)
    {
        //! if an external transformer is provided
        if (use_internal_transformer_)
        {
            delete sbt_;
            use_internal_transformer_ = false;
        }
        sbt_ = sbt;
    }
    else
    {
        // if no external transformer is provided
        if (!use_internal_transformer_)
        {
            sbt_ = new ModuleBase::SphericalBesselTransformer;
            use_internal_transformer_ = true;
        }
        // do nothing if an internal one is already in use
    }

    switch (update)
    {
    case 1:
        transform(true);
        break;
    case -1:
        transform(false);
        break;
    default:;
    }
}

void NumericalRadial::set_grid(const bool for_r_space, const int ngrid, const double* const grid, const char mode)
{
    assert(mode == 'i' || mode == 't');
    assert(ngrid > 1);

    double*& grid_to_set = (for_r_space ? rgrid_ : kgrid_);
    double*& value_to_set = (for_r_space ? rvalue_ : kvalue_);
    int& ngrid_to_set = (for_r_space ? nr_ : nk_);

    if (mode == 't')
    {

        // make sure a transform from the other space is available
        assert(for_r_space ? (kgrid_ && kvalue_) : (rgrid_ && rvalue_));

        delete[] grid_to_set;
        delete[] value_to_set;
        grid_to_set = new double[ngrid];
        value_to_set = new double[ngrid];
        ngrid_to_set = ngrid;

        for (int i = 0; i != ngrid; ++i)
        {
            grid_to_set[i] = grid[i];
        }

        check_fft_compliancy();
        transform(!for_r_space); // transform(true): r -> k; transform(false): k -> r
    }
    else
    {

        // cubic spline interpolation
        // step-1: compute f''(x[i])
        double* d2y = new double[ngrid_to_set];

        // FIXME boundary condition may not be correctly handled here
        // TODO should support more flexible boundary conditions or "not-a-knot" condition
        Mathzone_Add1::SplineD2(grid_to_set, value_to_set, ngrid_to_set, 0.0, 0.0, d2y);

        double* grid_tmp = new double[ngrid];
        double* value_tmp = new double[ngrid];
        double* dvalue_tmp = new double[ngrid];
        for (int i = 0; i != ngrid; ++i)
        {
            grid_tmp[i] = grid[i];
        }

        // step-2: interpolates values
        Mathzone_Add1::Cubic_Spline_Interpolation(grid_to_set,
                                                  value_to_set,
                                                  d2y,
                                                  ngrid_to_set,
                                                  grid_tmp,
                                                  ngrid,
                                                  value_tmp,
                                                  dvalue_tmp);
        delete[] d2y;
        delete[] dvalue_tmp;

        delete[] grid_to_set;
        delete[] value_to_set;
        grid_to_set = grid_tmp;
        value_to_set = value_tmp;
        ngrid_to_set = ngrid;

        check_fft_compliancy();
        transform(for_r_space); // transform(true): r -> k; transform(false): k -> r
    }
}

void NumericalRadial::set_uniform_grid(const bool for_r_space,
                                       const int ngrid,
                                       const double cutoff,
                                       const char mode,
                                       const bool enable_fft)
{
    double* grid = new double[ngrid];
    double dx = cutoff / (ngrid - 1);
    for (int i = 0; i != ngrid; ++i)
    {
        grid[i] = i * dx;
    }

    set_grid(for_r_space, ngrid, grid, mode);

    if (enable_fft)
    {
        set_uniform_grid(!for_r_space, ngrid, PI / dx, 't', false);
        is_fft_compliant_ = true;
    }

    delete[] grid;
}

void NumericalRadial::set_value(const bool for_r_space, const double* const value, const int p)
{
    if (for_r_space)
    {
        assert(rvalue_);
        for (int ir = 0; ir != nr_; ++ir)
        {
            rvalue_[ir] = value[ir];
        }
        pr_ = p;
        transform(true);
    }
    else
    {
        assert(kvalue_);
        for (int ik = 0; ik != nk_; ++ik)
        {
            kvalue_[ik] = value[ik];
        }
        pk_ = p;
        transform(false);
    }
}

void NumericalRadial::wipe(const bool r_space)
{
    if (r_space)
    {
        delete[] rgrid_;
        delete[] rvalue_;
        rgrid_ = nullptr;
        rvalue_ = nullptr;
        nr_ = 0;
        pr_ = 0;
    }
    else
    {
        delete[] kgrid_;
        delete[] kvalue_;
        kgrid_ = nullptr;
        kvalue_ = nullptr;
        nk_ = 0;
        pk_ = 0;
    }
    is_fft_compliant_ = false;
}

// TODO file format to be determined
// void NumericalRadial::save(std::string file) const {
//    if (file.empty()) {
//        file = symbol_ + "-" + std::to_string(l_) + "-" + std::to_string(ichi_) + ".dat";
//    }
//    std::ofstream ofs(file);
//    if (ofs.is_open()) {
//        ofs << "symbol " << symbol_ << std::endl;
//        ofs << "l " << l_ << std::endl;
//        ofs << "itype " << itype_ << std::endl;
//        ofs << "ichi " << ichi_ << std::endl;
//        ofs << "nr " << nr_ << std::endl;
//        ofs << "nk " << nk_ << std::endl;
//        ofs << "pr " << pr_ << std::endl;
//        ofs << "pk " << pk_ << std::endl;
//        ofs << "is_fft_compliant " << is_fft_compliant_ << std::endl;
//
//        if (rgrid_ && rvalue_) {
//            ofs << "rgrid " << "rvalue ";
//            for (int ir = 0; ir != nr_; ++ir) {
//                ofs << std::setw(20) << std::setprecision(15)
//                    << rgrid_[ir] << " "
//                    << rvalue_[ir] << std::endl;
//            }
//        }
//
//        ofs << std::endl;
//
//        if (kgrid_ && kvalue_) {
//            ofs << "kgrid " << "kvalue ";
//            for (int ik = 0; ik != nk_; ++ik) {
//                ofs << std::setw(20) << std::setprecision(15)
//                    << kgrid_[ik] << " " << kvalue_[ik] << std::endl;
//            }
//        }
//    }
//
//    ofs.close();
//
//}

void NumericalRadial::radtab(const char op,
                             const NumericalRadial& ket,
                             const int l,
                             double* const table,
                             const bool deriv)
{
    assert(op == 'S' || op == 'I' || op == 'T' || op == 'U');
    assert(l >= 0);

    // currently only FFT-compliant grids are supported!
    // FFT-based transform requires that two NumericalRadial objects have exactly the same grid
    assert(this->is_fft_compliant_ && ket.is_fft_compliant_);
    assert(this->nr_ == ket.nr_);
    assert(this->rcut() == ket.rcut());

    double* ktmp = new double[nk_];
    for (int ik = 0; ik != nk_; ++ik)
    {
        ktmp[ik] = this->kvalue_[ik] * ket.kvalue_[ik];
    }

    int op_pk = 0;
    switch (op)
    {
    case 'T':
        op_pk = -2;
        break;
    case 'U':
        op_pk = 2;
        break;
    default:;
    }

    if (deriv)
    { // derivative of the radial table
        if (l == 0)
        {
            // j'_0(x) = -j_1(x)
            sbt_->radrfft(1, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk - 1);
            for (int ir = 0; ir != nr_; ++ir)
            {
                table[ir] *= -1;
            }
        }
        else
        {
            // (2*l+1) * j'_l(x) = l * j_{l-1}(x) - (l+1) * j_{l+1}(x)
            double* rtmp = new double[nr_];
            sbt_->radrfft(l + 1, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk - 1);
            sbt_->radrfft(l - 1, nk_, this->kcut(), ktmp, rtmp, this->pk_ + ket.pk_ + op_pk - 1);
            for (int ir = 0; ir != nr_; ++ir)
            {
                table[ir] = (l * rtmp[ir] - (l + 1) * table[ir]) / (2 * l + 1);
            }
            delete[] rtmp;
        }
    }
    else
    {
        sbt_->radrfft(l, nk_, this->kcut(), ktmp, table, this->pk_ + ket.pk_ + op_pk);
    }
}

void NumericalRadial::transform(const bool forward)
{
    assert(forward ? (rgrid_ && rvalue_) : (kgrid_ && kvalue_));
    if ((forward && !kgrid_) || (!forward && !rgrid_))
    {
        return;
    }

    // currently only FFT-compliant grid is supported!
    assert(is_fft_compliant_);

    // value array must be pre-allocated
    if (forward)
    {
        sbt_->radrfft(l_, nr_, rgrid_[nr_ - 1], rvalue_, kvalue_, pr_);
        pk_ = 0;
    }
    else
    {
        sbt_->radrfft(l_, nk_, kgrid_[nk_ - 1], kvalue_, rvalue_, pk_);
        pr_ = 0;
    }
}

void NumericalRadial::check_fft_compliancy()
{
    is_fft_compliant_ = false;

    if (!rgrid_ || !kgrid_ || nr_ != nk_ || nr_ < 2)
    {
        return;
    }

    double tol = 4.0 * std::numeric_limits<double>::epsilon();

    double dr = rgrid_[nr_ - 1] / (nr_ - 1);
    double dk = kgrid_[nk_ - 1] / (nk_ - 1);
    if (std::abs(dr * dk - PI / (nr_ - 1)) > tol)
    {
        return;
    }

    for (int i = 0; i != nr_; ++i)
    {
        if (std::abs(i * dr - rgrid_[i]) > tol || std::abs(i * dk - kgrid_[i]) > tol)
        {
            return;
        }
    }

    is_fft_compliant_ = true;
}
