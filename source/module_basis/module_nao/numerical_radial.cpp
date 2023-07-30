#include "module_basis/module_nao/numerical_radial.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <numeric>

#include "module_base/constants.h"
#include "module_base/cubic_spline.h"
#include "module_base/math_integral.h"
#include "module_base/spherical_bessel_transformer.h"

using ModuleBase::PI;

NumericalRadial::NumericalRadial() :
    sbt_(new ModuleBase::SphericalBesselTransformer),
    use_internal_transformer_(true)
{
}

NumericalRadial::NumericalRadial(const NumericalRadial& other) :
    symbol_(other.symbol_),
    itype_(other.itype_),
    l_(other.l_),
    izeta_(other.izeta_),
    nr_(other.nr_),
    nk_(other.nk_),
    is_fft_compliant_(other.is_fft_compliant_),
    pr_(other.pr_),
    pk_(other.pk_),
    sbt_(other.use_internal_transformer_ ? new ModuleBase::SphericalBesselTransformer : other.sbt_),
    use_internal_transformer_(other.use_internal_transformer_)
{
    // deep copy
    if (other.ptr_rgrid())
    {
        rgrid_ = new double[nr_];
        rvalue_ = new double[nr_];
        std::memcpy(rgrid_, other.rgrid_, nr_ * sizeof(double));
        std::memcpy(rvalue_, other.rvalue_, nr_ * sizeof(double));
    }

    if (other.ptr_kgrid())
    {
        kgrid_ = new double[nk_];
        kvalue_ = new double[nk_];
        std::memcpy(kgrid_, other.kgrid_, nk_ * sizeof(double));
        std::memcpy(kvalue_, other.kvalue_, nk_ * sizeof(double));
    }
}

NumericalRadial& NumericalRadial::operator=(const NumericalRadial& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    // wipe off r & k space data
    wipe(true);
    wipe(false);

    symbol_ = rhs.symbol_;
    itype_ = rhs.itype_;
    izeta_ = rhs.izeta_;
    l_ = rhs.l_;

    nr_ = rhs.nr_;
    nk_ = rhs.nk_;

    is_fft_compliant_ = rhs.is_fft_compliant_;

    pr_ = rhs.pr_;
    pk_ = rhs.pk_;

    // deep copy
    if (rhs.ptr_rgrid())
    {
        rgrid_ = new double[nr_];
        rvalue_ = new double[nr_];
        std::memcpy(rgrid_, rhs.rgrid_, nr_ * sizeof(double));
        std::memcpy(rvalue_, rhs.rvalue_, nr_ * sizeof(double));
    }

    if (rhs.ptr_kgrid())
    {
        kgrid_ = new double[nk_];
        kvalue_ = new double[nk_];
        std::memcpy(kgrid_, rhs.kgrid_, nk_ * sizeof(double));
        std::memcpy(kvalue_, rhs.kvalue_, nk_ * sizeof(double));
    }

    set_transformer(rhs.use_internal_transformer_ ? nullptr : rhs.sbt_, 0);

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
                            const int izeta,
                            const std::string symbol,
                            const int itype)
{
    assert(l >= 0);
    assert(ngrid > 1);
    assert(grid && value);

    // grid must be strictly increasing and every element must be non-negative
    assert(std::is_sorted(grid, grid + ngrid, std::less_equal<double>())); // std::less<>() would allow equal values
    assert(grid[0] >= 0.0);

    // wipe off any existing r & k space data
    wipe(true);
    wipe(false);

    symbol_ = symbol;
    itype_ = itype;
    izeta_ = izeta;
    l_ = l;

    if (for_r_space)
    {
        nr_ = ngrid;
        pr_ = p;
        rgrid_ = new double[nr_];
        rvalue_ = new double[nr_];
        std::memcpy(rgrid_, grid, nr_ * sizeof(double));
        std::memcpy(rvalue_, value, nr_ * sizeof(double));
    }
    else
    {
        nk_ = ngrid;
        pk_ = p;
        kgrid_ = new double[nk_];
        kvalue_ = new double[nk_];
        std::memcpy(kgrid_, grid, nk_ * sizeof(double));
        std::memcpy(kvalue_, value, nk_ * sizeof(double));
    }
}

void NumericalRadial::set_transformer(ModuleBase::SphericalBesselTransformer* sbt, int update)
{
    assert(update == 0 || update == 1 || update == -1);

    if (use_internal_transformer_ && sbt)
    { // internal -> external
        delete sbt_;
        use_internal_transformer_ = false;
        sbt_ = sbt;
    }
    else if (!use_internal_transformer_ && !sbt)
    { // external -> internal
        sbt_ = new ModuleBase::SphericalBesselTransformer;
        use_internal_transformer_ = true;
    }
    else if (!use_internal_transformer_ && sbt)
    { // external -> another external
        sbt_ = sbt;
    }

    switch (update)
    {
    case 1:
        transform(true); // forward transform r -> k
        break;
    case -1:
        transform(false); // backward transform k -> r
        break;
    default:;
    }
}

void NumericalRadial::set_grid(const bool for_r_space, const int ngrid, const double* const grid, const char mode)
{
    assert(mode == 'i' || mode == 't');
    assert(ngrid > 1);

    // grid must be strictly increasing and every element must be non-negative
    assert(std::is_sorted(grid, grid + ngrid, std::less_equal<double>())); // std::less<>() would allow equal values
    assert(grid[0] >= 0.0);

    // tbu stands for "to be updated"
    double*& grid_tbu = (for_r_space ? rgrid_ : kgrid_);
    double*& value_tbu = (for_r_space ? rvalue_ : kvalue_);
    int& ngrid_tbu = (for_r_space ? nr_ : nk_);

    if (mode == 't')
    { // obtain new values by a transform from the other space
        // make sure a transform from the other space is available
        assert(for_r_space ? (kgrid_ && kvalue_) : (rgrid_ && rvalue_));

        delete[] grid_tbu;
        delete[] value_tbu;
        grid_tbu = new double[ngrid];
        value_tbu = new double[ngrid];
        ngrid_tbu = ngrid;
        std::memcpy(grid_tbu, grid, ngrid * sizeof(double));

        is_fft_compliant_ = is_fft_compliant(nr_, rgrid_, nk_, kgrid_);
        transform(!for_r_space); // transform(true): r -> k; transform(false): k -> r
    }
    else
    { // obtain new values by interpolation in the current space
        // make sure an interpolation in the current space is available
        assert(grid_tbu && value_tbu);

        // cubic spline interpolation
        ModuleBase::CubicSpline cubspl;
        cubspl.build(ngrid_tbu, grid_tbu, value_tbu); // not-a-knot boundary condition

        double* grid_new = new double[ngrid];
        double* value_new = new double[ngrid];

        std::memcpy(grid_new, grid, ngrid * sizeof(double));
        std::fill_n(value_new, ngrid, 0.0);

        // do interpolation for grid points within the range of the origional grid
        // for grid points outside the original range, simply set the values to zero

        // grid_start is the first grid point that is greater than or equal to grid_tbu[0]
        double* grid_start = std::lower_bound(grid_new, grid_new + ngrid, grid_tbu[0]);

        // grid_end is the first grid point that is strictly greater than grid_tbu[ngrid_tbu-1]
        double* grid_end = std::upper_bound(grid_new, grid_new + ngrid, grid_tbu[ngrid_tbu - 1]);

        cubspl.interp(std::distance(grid_start, grid_end), grid_start, value_new + std::distance(grid_new, grid_start));

        delete[] grid_tbu;
        delete[] value_tbu;

        grid_tbu = grid_new;
        value_tbu = value_new;
        ngrid_tbu = ngrid;

        is_fft_compliant_ = is_fft_compliant(nr_, rgrid_, nk_, kgrid_);
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
    delete[] grid;

    if (enable_fft)
    {
        set_uniform_grid(!for_r_space, ngrid, PI / dx, 't', false);
    }
}

void NumericalRadial::set_value(const bool for_r_space, const double* const value, const int p)
{
    if (for_r_space)
    {
        assert(rvalue_);
        std::memcpy(rvalue_, value, nr_ * sizeof(double));
        pr_ = p;
        transform(true);
    }
    else
    {
        assert(kvalue_);
        std::memcpy(kvalue_, value, nk_ * sizeof(double));
        pk_ = p;
        transform(false);
    }
}

void NumericalRadial::wipe(const bool r_space)
{
    // wipe the grid and value in the given space
    // if r_space is true, wipe r space; otherwise wipe k space
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
//        file = symbol_ + "-" + std::to_string(l_) + "-" + std::to_string(izeta_) + ".dat";
//    }
//    std::ofstream ofs(file);
//    if (ofs.is_open()) {
//        ofs << "symbol " << symbol_ << std::endl;
//        ofs << "l " << l_ << std::endl;
//        ofs << "itype " << itype_ << std::endl;
//        ofs << "izeta " << izeta_ << std::endl;
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
                             const int nr_tab,
                             const double* const rgrid_tab,
                             const bool deriv) const
{
    assert(op == 'S' || op == 'I' || op == 'T' || op == 'U');
    assert(l >= 0);
    assert(rgrid_tab && nr_tab > 0);

    // radtab requires that two NumericalRadial objects have exactly the same (non-null) kgrid_
    assert(nk_ > 0 && nk_ == ket.nk_);
    assert(std::equal(kgrid_, kgrid_ + nk_, ket.kgrid_));

    bool use_radrfft = is_fft_compliant(nr_tab, rgrid_tab, nk_, kgrid_);

    // function to undergo a spherical Bessel transform:
    // overlap: chi1(k) * chi2(k)
    // kinetic: k^2 * chi1(k) * chi2(k)
    // Coulomb: k^(-2) * chi1(k) * chi2(k)
    double* fk = new double[nk_];
    std::transform(kvalue_, kvalue_ + nk_, ket.kvalue_, fk, std::multiplies<double>());

    int op_pk = 0;
    switch (op)
    {
    case 'T':
        op_pk = -2;
        break;
    case 'U':
        op_pk = 2;
        break;
    default:; // for overlap integral op_pk = 0
    }

    if (use_radrfft)
    {
        sbt_->radrfft(l, nk_, kcut(), fk, table, pk_ + ket.pk_ + op_pk, deriv);
    }
    else
    {
        sbt_->direct(l, nk_, kgrid_, fk, nr_tab, rgrid_tab, table, pk_ + ket.pk_ + op_pk, deriv);
    }

    delete[] fk;

    // spherical Bessel transform has a prefactor of sqrt(2/pi) while the prefactor for the radial table
    // of two-center integrals is 4*pi
    double pref = ModuleBase::FOUR_PI * std::sqrt(ModuleBase::PI / 2.0);
    std::for_each(table, table + nr_tab, [pref](double& x) { x *= pref; });
}

void NumericalRadial::normalize(bool for_r_space)
{
    int& ngrid = for_r_space ? nr_ : nk_;

    // tbu stands for "to be updated"
    double*& grid_tbu = for_r_space ? rgrid_ : kgrid_;
    double*& value_tbu = for_r_space ? rvalue_ : kvalue_;

    double factor = 0.0;
    double* integrand = new double[ngrid];
    double* rab = new double[ngrid];

    std::adjacent_difference(grid_tbu, grid_tbu + ngrid, rab);
    std::transform(value_tbu, value_tbu + ngrid, grid_tbu, integrand, std::multiplies<double>());
    std::for_each(integrand, integrand + ngrid, [](double& x) { x *= x; });

    factor = ModuleBase::Integral::simpson(ngrid, integrand, &rab[1]);
    factor = 1. / std::sqrt(factor);

    std::for_each(value_tbu, value_tbu + ngrid, [factor](double& x) { x *= factor; });
    transform(for_r_space);
    delete[] rab;
    delete[] integrand;
}

void NumericalRadial::transform(const bool forward)
{
    // grid & value must exist in the initial space
    assert(forward ? (rgrid_ && rvalue_) : (kgrid_ && kvalue_));

    // do nothing if there is no grid in the destination space
    if ((forward && !kgrid_) || (!forward && !rgrid_))
    {
        return;
    }

    if (forward)
    { // r -> k
        if (is_fft_compliant_)
        {
            sbt_->radrfft(l_, nr_, rgrid_[nr_ - 1], rvalue_, kvalue_, pr_);
        }
        else
        {
            sbt_->direct(l_, nr_, rgrid_, rvalue_, nk_, kgrid_, kvalue_, pr_);
        }
        pk_ = 0;
    }
    else
    { // k -> r
        if (is_fft_compliant_)
        {
            sbt_->radrfft(l_, nk_, kgrid_[nk_ - 1], kvalue_, rvalue_, pk_);
        }
        else
        {
            sbt_->direct(l_, nk_, kgrid_, kvalue_, nr_, rgrid_, rvalue_, pk_);
        }
        pr_ = 0;
    }
}

bool NumericalRadial::is_fft_compliant(const int nr,
                                       const double* const rgrid,
                                       const int nk,
                                       const double* const kgrid) const
{
    if (!rgrid || !kgrid || nr != nk || nr < 2)
    {
        return false;
    }

    double dr = rgrid[nr - 1] / (nr - 1);
    double dk = kgrid[nk - 1] / (nk - 1);
    double tol = 4.0 * std::numeric_limits<double>::epsilon();

    return std::abs(dr * dk - PI / (nr - 1)) < tol
           && std::all_of(rgrid, rgrid + nr, [&](const double& r) { return std::abs(r - (&r - rgrid) * dr) < tol; })
           && std::all_of(kgrid, kgrid + nk, [&](const double& k) { return std::abs(k - (&k - kgrid) * dk) < tol; });
}
