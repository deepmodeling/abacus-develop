#include "module_basis/module_nao/two_center_integral.h"

#include "module_base/ylm.h"

TwoCenterIntegral::TwoCenterIntegral() : with_deriv_(false), use_internal_gaunt_(false), rgt_(nullptr)
{
}

TwoCenterIntegral::~TwoCenterIntegral()
{
    if (use_internal_gaunt_)
    {
        delete rgt_;
    }
}

void TwoCenterIntegral::build(const RadialCollection& bra,
                              const RadialCollection& ket,
                              const char op,
                              const int nr,
                              const double cutoff,
                              const bool with_deriv,
                              RealGauntTable* const rgt)
{
    with_deriv_ = with_deriv;
    table_.build(bra, ket, op, nr, cutoff, with_deriv);

    if (rgt)
    { // if an external gaunt table is provided
        if (use_internal_gaunt_)
        {
            delete rgt_;    
        }
        rgt_ = rgt;
    }
    else
    { // if no external gaunt table is provided (which implies an internal one)
        if (!use_internal_gaunt_)
        {
            rgt_ = new RealGauntTable;
        }
    }

    rgt_->build(std::max(bra.lmax(), ket.lmax()));
}

void TwoCenterIntegral::get(const int itype1, 
                            const int l1, 
                            const int izeta1, 
                            const int m1, 
                            const int itype2,
                            const int l2,
                            const int izeta2,
                            const int m2,
	                        const ModuleBase::Vector3<double>& vR, // R = R2 - R1
                            const bool deriv,
                            double* out) const
{
    assert( (deriv && with_deriv_) || !deriv );

    double R = vR.norm();
    std::fill(out, out + (deriv ? 3 : 1), 0.0);

    // generate all necessary spherical harmonics
	std::vector<double> Rl_Y;
	std::vector<std::vector<double>> grad_Rl_Y;

    int lmax = l1 + l2;
	if (deriv)
	{
		ModuleBase::Ylm::grad_rl_sph_harm(lmax, vR[0], vR[1], vR[2], Rl_Y, grad_Rl_Y);
	}
	else
	{
		ModuleBase::Ylm::rl_sph_harm(lmax, vR[0], vR[1], vR[2], Rl_Y);
	}

    int sign = 1;
    double tiny = 1e-8;
    for (int l = std::abs(l1 - l2); l <= l1 + l2; l += 2)
    {
        double S_by_Rl = table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, false);
        double tmp; // ( (dS/dR)/R^l - l/R * (S/R^l)  ) / R
        if (deriv)
        {
            tmp = R < tiny ? 0.0 : 
                             ( table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, true) / std::pow(R, l) - l / R * S_by_Rl ) / R;
        }

		for (int m = -l; m < l; ++m)
        {
            double G = (*rgt_)(l1, l2, l, m1, m2, m);

            if (deriv)
            {
                for (int i = 0; i < 3; ++i)
                {
                    out[i] += sign * G * (tmp * vR[i] * Rl_Y[lm_index(l, m)] + S_by_Rl * grad_Rl_Y[lm_index(l, m)][i]);
                }
            }
            else
            {
                out[0] += sign * G * S_by_Rl * Rl_Y[lm_index(l, m)];
            }

        }
        sign = -sign;
    }
}

int TwoCenterIntegral::lm_index(const int l, const int m) const
{
    return l * l + (m > 0 ? 2 * m - 1 : -2 * m);
}
