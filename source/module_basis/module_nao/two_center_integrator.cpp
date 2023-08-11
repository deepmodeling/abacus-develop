#include "module_basis/module_nao/two_center_integrator.h"

#include "module_base/vector3.h"
#include "module_base/ylm.h"

TwoCenterIntegrator::TwoCenterIntegrator():
    is_tabulated_(false),
    op_('\0'),
    with_deriv_(false)
{
}

void TwoCenterIntegrator::tabulate(const RadialCollection& bra,
                                 const RadialCollection& ket,
                                 const char op,
                                 const int nr,
                                 const double cutoff,
                                 const bool with_deriv)
{
    op_ = op;
    with_deriv_ = with_deriv;
    table_.build(bra, ket, op, nr, cutoff, with_deriv);
    RealGauntTable::instance().build(std::max(bra.lmax(), ket.lmax()));
    is_tabulated_ = true;
}

void TwoCenterIntegrator::calculate(const int itype1, 
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
    assert( is_tabulated_ );
    assert( (deriv && with_deriv_) || !deriv );

    double R = vR.norm();

    // unit vector along R
    ModuleBase::Vector3<double> uR = (R == 0.0 ? ModuleBase::Vector3<double>(0., 0., 1.) : vR / R);

    std::fill(out, out + (deriv ? 3 : 1), 0.0);

    // generate all necessary real spherical harmonics (multiplied by R^l)
	std::vector<double> Rl_Y;
	std::vector<std::vector<double>> grad_Rl_Y;

	if (deriv)
	{
		ModuleBase::Ylm::grad_rl_sph_harm(l1 + l2, vR[0], vR[1], vR[2], Rl_Y, grad_Rl_Y);
	}
	else
	{
		ModuleBase::Ylm::rl_sph_harm(l1 + l2, vR[0], vR[1], vR[2], Rl_Y);
	}

    int sign = 1;
    for (int l = std::abs(l1 - l2); l <= l1 + l2; l += 2)
    {
        double S_by_Rl = table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, false);

        // tmp is (d/dR)(S/R^l)
        double tmp = (deriv && R > 1e-6) ? 
            table_.lookup(itype1, l1, izeta1, itype2, l2, izeta2, l, R, true) / std::pow(R, l)
                - l / R * S_by_Rl
            : 0.0;

		for (int m = -l; m <= l; ++m)
        {
            double G = RealGauntTable::instance()(l1, l2, l, m1, m2, m);

            if (deriv)
            {
                for (int i = 0; i < 3; ++i)
                {
                    out[i] += sign * G * ( tmp * uR[i] * Rl_Y[ylm_index(l, m)] 
                                            + S_by_Rl * grad_Rl_Y[ylm_index(l, m)][i] );
                }
            }
            else
            {
                out[0] += sign * G * S_by_Rl * Rl_Y[ylm_index(l, m)];
            }
        }
        sign = -sign;
    }
}

int TwoCenterIntegrator::ylm_index(const int l, const int m) const
{
    return l * l + (m > 0 ? 2 * m - 1 : -2 * m);
}
