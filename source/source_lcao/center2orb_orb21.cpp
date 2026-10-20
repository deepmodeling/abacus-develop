//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#include "center2orb_orb21.h"

#include "source_base/constants.h"
#include "source_base/ylm.h"

#include <cmath>

Center2_Orb::Orb21::Orb21(const Numerical_Orbital_Lm& nA1_in,
                          const Numerical_Orbital_Lm& nA2_in,
                          const Numerical_Orbital_Lm& nB_in,
                          const ModuleBase::Sph_Bessel_Recursive::D2* psb,
                          const ORB_gaunt_table& MGT_in)
    : nA1(nA1_in), nA2(nA2_in), nB(nB_in), psb_(psb), MGT(MGT_in)
{
}

void Center2_Orb::Orb21::prepare_radial_table()
{
    auto products = std::make_shared<ProductOrbitals>();
    const Numerical_Orbital_Lm& nA_short = (this->nA1.getNr() <= this->nA2.getNr()) ? this->nA1 : this->nA2;

    std::vector<double> nA_tmp(nA_short.getNr());
    for (size_t ir = 0; ir != nA_tmp.size(); ++ir)
    {
        nA_tmp[ir] = this->nA1.getPsi(ir) * this->nA2.getPsi(ir);
    }

    const int LA1 = this->nA1.getL();
    const int LA2 = this->nA2.getL();
    for (int LA = std::abs(LA1 - LA2); LA <= LA1 + LA2; ++LA)
    {
        if ((LA - std::abs(LA1 - LA2)) % 2 == 1) // if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
        {
            continue;
        }

        (*products)[LA].set_orbital_info(nA_short.getLabel(),
                                         nA_short.getType(),
                                         LA,
                                         1, // N?
                                         nA_short.getNr(),
                                         nA_short.getRab(),
                                         nA_short.getRadial(),
                                         Numerical_Orbital_Lm::Psi_Type::Psi,
                                         nA_tmp.data(),
                                         nA_short.getNk(),
                                         nA_short.getDk(),
                                         nA_short.getDruniform(),
                                         false,
                                         true,
                                         false); // Product orbitals use reciprocal-space gradients, not zty.
    }
    this->nA = std::move(products);
    this->bind_product_orbitals();
}

void Center2_Orb::Orb21::bind_product_orbitals()
{
    this->orb11s.clear();
    for (const auto& product: *this->nA)
    {
        this->orb11s.emplace(product.first, Orb11(product.second, this->nB, this->psb_, this->MGT));
    }
}

void Center2_Orb::Orb21::init_radial_table()
{
    this->prepare_radial_table();
    for (auto& orbital: this->orb11s)
    {
        orbital.second.init_radial_table();
    }
}

void Center2_Orb::Orb21::init_radial_table(const std::set<size_t>& radials)
{
    this->prepare_radial_table();
    for (auto& orbital: this->orb11s)
    {
        orbital.second.init_radial_table(radials);
    }
}

void Center2_Orb::Orb21::init_radial_table(const std::vector<Orb21*>& orbitals, const bool use_gpu)
{
    const auto pairs = prepare_radial_tables(orbitals);
    Orb11::init_radial_table(pairs, use_gpu);
}

void Center2_Orb::Orb21::init_radial_table(const std::vector<Orb21*>& orbitals, const std::set<size_t>& radials)
{
    const auto pairs = prepare_radial_tables(orbitals);
    for (Orb11* pair: pairs)
    {
        pair->init_radial_table(radials);
    }
}

std::vector<Center2_Orb::Orb11*> Center2_Orb::Orb21::prepare_radial_tables(const std::vector<Orb21*>& orbitals)
{
    std::vector<Orb11*> pairs;
    // A product depends on its two orbitals on center A, not the orbital on
    // center B. Share it within this build and retain ownership with each user.
    std::map<const Numerical_Orbital_Lm*, std::map<const Numerical_Orbital_Lm*, std::shared_ptr<const ProductOrbitals>>>
        products;
    for (Orb21* orbital: orbitals)
    {
        // Orbital transforms grow the shared Bessel cache, so prepare them serially.
        auto& product = products[&orbital->nA1][&orbital->nA2];
        if (!product)
        {
            orbital->prepare_radial_table();
            product = orbital->nA;
        }
        else
        {
            orbital->nA = product;
            orbital->bind_product_orbitals();
        }
        for (auto& pair: orbital->orb11s)
        {
            pairs.push_back(&pair.second);
        }
    }
    return pairs;
}

double Center2_Orb::Orb21::cal_overlap(const ModuleBase::Vector3<double>& RA,
                                       const ModuleBase::Vector3<double>& RB,
                                       const int& mA1,
                                       const int& mA2,
                                       const int& mB) const
{
    const int LA1 = this->nA1.getL();
    const int LA2 = this->nA2.getL();
    const int idx1 = this->MGT.get_lm_index(LA1, mA1);
    const int idx2 = this->MGT.get_lm_index(LA2, mA2);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    double overlap = 0.0;

    for (const auto& orb11: this->orb11s)
    {
        const int LA = orb11.first;

        for (int mA = 0; mA != 2 * LA + 1; ++mA)
        // const int mA=mA1+mA2;
        {
            const double Gaunt_real_A1_A2_A12 = *(Gaunt_Coefficients_ptr + this->MGT.get_lm_index(LA, mA));
            if (0 == Gaunt_real_A1_A2_A12)
            {
                continue;
            }

            overlap += Gaunt_real_A1_A2_A12 * orb11.second.cal_overlap(RA, RB, mA, mB);
        }
    }

    return overlap;
}

ModuleBase::Vector3<double> Center2_Orb::Orb21::cal_grad_overlap(const ModuleBase::Vector3<double>& RA,
                                                                 const ModuleBase::Vector3<double>& RB,
                                                                 const int& mA1,
                                                                 const int& mA2,
                                                                 const int& mB) const
{
    const int LA1 = this->nA1.getL();
    const int LA2 = this->nA2.getL();
    const int idx1 = this->MGT.get_lm_index(LA1, mA1);
    const int idx2 = this->MGT.get_lm_index(LA2, mA2);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    ModuleBase::Vector3<double> grad_overlap(0.0, 0.0, 0.0);

    for (const auto& orb11: this->orb11s)
    {
        const int LA = orb11.first;

        for (int mA = 0; mA != 2 * LA + 1; ++mA)
        // const int mA=mA1+mA2;
        {
            const double Gaunt_real_A1_A2_A12 = *(Gaunt_Coefficients_ptr + this->MGT.get_lm_index(LA, mA));
            if (0 == Gaunt_real_A1_A2_A12)
            {
                continue;
            }

            grad_overlap += Gaunt_real_A1_A2_A12 * orb11.second.cal_grad_overlap(RA, RB, mA, mB);
        }
    }

    return grad_overlap;
}

template <typename T>
std::vector<T> Center2_Orb::Orb21::cal_magnetic_components(
    const ModuleBase::Vector3<double>& RA,
    const ModuleBase::Vector3<double>& RB,
    std::vector<T> (Orb11::*evaluate)(const ModuleBase::Vector3<double>&,
                                    const ModuleBase::Vector3<double>&) const) const
{
    const int LA1 = this->nA1.getL();
    const int LA2 = this->nA2.getL();
    const int count_A1 = 2 * LA1 + 1;
    const int count_A2 = 2 * LA2 + 1;
    const int count_B = 2 * this->nB.getL() + 1;
    const int count = count_A1 * count_A2 * count_B;
    std::vector<T> result(count);
    for (const auto& orbital: this->orb11s)
    {
        const int LA = orbital.first;
        const int count_A = 2 * LA + 1;
        // These two-center overlaps are shared by all ways to couple the
        // magnetic components of the two original orbitals on center A.
        const std::vector<T> overlaps = (orbital.second.*evaluate)(RA, RB);
        for (int ma1 = 0; ma1 < count_A1; ++ma1)
        {
            const int idx1 = this->MGT.get_lm_index(LA1, ma1);
            for (int ma2 = 0; ma2 < count_A2; ++ma2)
            {
                const int idx2 = this->MGT.get_lm_index(LA2, ma2);
                const int offset = (ma1 * count_A2 + ma2) * count_B;
                for (int ma = 0; ma < count_A; ++ma)
                {
                    const int idx = this->MGT.get_lm_index(LA, ma);
                    const double coefficient = this->MGT.Gaunt_Coefficients(idx1, idx2, idx);
                    if (coefficient == 0.0)
                    {
                        continue;
                    }
                    for (int mb = 0; mb < count_B; ++mb)
                    {
                        result[offset + mb] += coefficient * overlaps[ma * count_B + mb];
                    }
                }
            }
        }
    }
    return result;
}

std::vector<double> Center2_Orb::Orb21::cal_overlap(const ModuleBase::Vector3<double>& RA,
                                                    const ModuleBase::Vector3<double>& RB) const
{
    return this->cal_magnetic_components<double>(RA, RB, &Orb11::cal_overlap);
}

std::vector<ModuleBase::Vector3<double>> Center2_Orb::Orb21::cal_grad_overlap(
    const ModuleBase::Vector3<double>& RA,
    const ModuleBase::Vector3<double>& RB) const
{
    return this->cal_magnetic_components<ModuleBase::Vector3<double>>(RA, RB, &Orb11::cal_grad_overlap);
}
