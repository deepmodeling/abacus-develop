#include "center2_orb.h"

#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/mathzone_add1.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

int Center2_Orb::get_rmesh(const double& R1, const double& R2, const double dr)
{
    int rmesh = static_cast<int>((R1 + R2) / dr) + 5;
    // mohan update 2009-09-08 +1 ==> +5
    // considering interpolation or so on...
    if (rmesh % 2 == 0) {
        rmesh++;
    }

    if (rmesh <= 0)
    {
        // GlobalV::ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
        // GlobalV::ofs_warning << "\n rmesh = " << rmesh;
        std::cout << "\n R1 = " << R1 << " R2 = " << R2;
        std::cout << "\n rmesh = " << rmesh;
        ModuleBase::WARNING_QUIT("Center2_Orb::get_rmesh", "rmesh <= 0");
    }
    return rmesh;
}

// Peize Lin update 2016-01-26
void Center2_Orb::init_Table_Spherical_Bessel(const int Lmax_used,
                                              const double dr,
                                              const double dk,
                                              const int kmesh,
                                              const int Rmesh,
                                              ModuleBase::Sph_Bessel_Recursive::D2*& psb)
{
    ModuleBase::TITLE("Center2_Orb", "init_Table_Spherical_Bessel");

    for (auto& sb: ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool)
    {
        if (dr * dk == sb.get_dx())
        {
            psb = &sb;
            break;
        }
    }

    if (!psb)
    {
        ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.push_back({});
        psb = &ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.back();
    }

    psb->set_dx(dr * dk);
    psb->cal_jlx(Lmax_used+1, Rmesh, kmesh);                // +1 for drs needs psb.jlx[l+1]. Peize Lin update 2025-12-27
}

// Peize Lin accelerate 2017-10-02
void Center2_Orb::cal_ST_Phi12_R(const int& job,
                                 const int& l,
                                 const Numerical_Orbital_Lm& n1,
                                 const Numerical_Orbital_Lm& n2,
                                 const int& rmesh,
                                 std::vector<double> &rs,
                                 std::vector<double> &drs,
                                 const ModuleBase::Sph_Bessel_Recursive::D2* psb)
{
    ModuleBase::timer::tick("Center2_Orb", "cal_ST_Phi12_R");

    assert(rmesh <= rs.size());
    assert(rmesh <= drs.size());

    const int kmesh = n1.getNk();
    const double* kpoint = n1.getKpoint();
    const double dk = n1.getDk();

    std::vector<double> k1_dot_k2(kmesh);
    // Peize Lin change 2017-12-12
    switch (job)
    {
    case 1: // calculate overlap
        if (!n1.get_psif().empty() && !n2.get_psi_k2().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsif(ik) * n2.getPsi_k2(ik);
            }
        }
        else if (!n1.get_psi_k().empty() && !n2.get_psi_k().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
            }
        }
        else if (!n1.get_psi_k2().empty() && !n2.get_psif().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsif(ik);
            }
        }
        break;

    case 2: // calculate kinetic energy
        for (int ik = 0; ik < kmesh; ik++)
        {
            k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsi_k2(ik);
        }
        break;
    }

    std::vector<double> k1_dot_k2_dot_kpoint(kmesh);
    for (int ik = 0; ik < kmesh; ik++)
    {
        k1_dot_k2_dot_kpoint[ik] = k1_dot_k2[ik] * kpoint[ik];
    }

    // Drs
    // djl = (l*j(l-1) - (l+1)j(l+1))/(2l+1)

    // previous version

    // double* integrated_func = new double[kmesh];

    assert(psb->get_jlx().size()>=l+2);
    const int lml = (l>0) ? (l-1) : 0;
    const std::vector<std::vector<double>>& jlm1 = psb->get_jlx().at(lml);
    const std::vector<std::vector<double>>& jl = psb->get_jlx().at(l);
    const std::vector<std::vector<double>>& jlp1 = psb->get_jlx().at(l+1);
    assert(jlm1.size()>=rmesh);
    assert(jl.size()>=rmesh);
    assert(jlp1.size()>=rmesh);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ir = 0; ir < rmesh; ir++)
    {
        std::vector<double> integrated_func(kmesh);
        const std::vector<double>& jl_r = jl[ir];
        assert(jl_r.size()>=kmesh);
        for (int ik = 0; ik < kmesh; ++ik)
        {
            integrated_func[ik] = jl_r[ik] * k1_dot_k2[ik];
        }
        // Call simpson integration
        double temp = 0.0;

        ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), dk, temp);
        rs[ir] = temp * ModuleBase::FOUR_PI;

        // Peize Lin accelerate 2017-10-02
        const std::vector<double>& jlm1_r = jlm1[ir];
        const std::vector<double>& jlp1_r = jlp1[ir];
        assert(jlm1_r.size()>=kmesh);
        assert(jlp1_r.size()>=kmesh);
        const double fac = l / (l + 1.0);
        if (l == 0)
        {
            for (int ik = 0; ik < kmesh; ++ik)
            {
                integrated_func[ik] = jlp1_r[ik] * k1_dot_k2_dot_kpoint[ik];
            }
        }
        else
        {
            for (int ik = 0; ik < kmesh; ++ik)
            {
                integrated_func[ik] = (jlp1_r[ik] - fac * jlm1_r[ik]) * k1_dot_k2_dot_kpoint[ik];
            }
        }

        ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), dk, temp);
        drs[ir] = -ModuleBase::FOUR_PI * (l + 1) / (2.0 * l + 1) * temp;
    }

    // liaochen modify on 2010/4/22
    // special case for R=0
    // we store Slm(R) / R**l at the fisrt point, rather than Slm(R)

    if (l > 0)
    {
        std::vector<double> integrated_func(kmesh);
        double temp = 0.0;

        for (int ik = 0; ik < kmesh; ik++)
        {
            integrated_func[ik] = k1_dot_k2[ik] * std::pow(kpoint[ik], l);
        }

        ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), dk, temp);
        rs[0] = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac(2 * l + 1) * temp;
    }

    ModuleBase::timer::tick("Center2_Orb", "cal_ST_Phi12_R");
}

#include "source_base/constants.h"

// Peize Lin add 2017-10-27
void Center2_Orb::cal_ST_Phi12_R(const int& job,
                                 const int& l,
                                 const Numerical_Orbital_Lm& n1,
                                 const Numerical_Orbital_Lm& n2,
                                 const std::set<size_t>& radials,
                                 std::vector<double> &rs,
                                 std::vector<double> &drs,
                                 const ModuleBase::Sph_Bessel_Recursive::D2* psb)
{
    //	ModuleBase::TITLE("Center2_Orb","cal_ST_Phi12_R");
    ModuleBase::timer::tick("Center2_Orb", "cal_ST_Phi12_R");

    const int kmesh = n1.getNk();
    const double* kpoint = n1.getKpoint();
    const double dk = n1.getDk();

    std::vector<double> k1_dot_k2(kmesh);
    switch (job)
    {
    case 1: // calculate overlap
        if (!n1.get_psif().empty() && !n2.get_psi_k2().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsif(ik) * n2.getPsi_k2(ik);
            }
        }
        else if (!n1.get_psi_k().empty() && !n2.get_psi_k().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
            }
        }
        else if (!n1.get_psi_k2().empty() && !n2.get_psif().empty())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsif(ik);
            }
        }
        break;

    case 2: // calculate kinetic energy
        for (int ik = 0; ik < kmesh; ik++)
        {
            k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsi_k2(ik);
        }
        break;
    }

    std::vector<double> k1_dot_k2_dot_kpoint(kmesh);
    for (int ik = 0; ik < kmesh; ik++)
    {
        k1_dot_k2_dot_kpoint[ik] = k1_dot_k2[ik] * kpoint[ik];
    }

    std::vector<double> integrated_func(kmesh);

    assert(psb->get_jlx().size()>=l+2);
    const int lm1 = (l>0) ? (l-1) : 0;
    const std::vector<std::vector<double>>& jlm1 = psb->get_jlx().at(lm1);
    const std::vector<std::vector<double>>& jl = psb->get_jlx().at(l);
    const std::vector<std::vector<double>>& jlp1 = psb->get_jlx().at(l+1);

    for (const size_t& ir: radials)
    {
        // if(rs[ir])  => rs[ir]  has been calculated
        // if(drs[ir]) => drs[ir] has been calculated
        // Actually, if(ir[ir]||dr[ir]) is enough. Double insurance for the sake of avoiding numerical errors
        if (rs.at(ir) && drs.at(ir)) {
            continue;
        }

        const std::vector<double>& jl_r = jl[ir];
        assert(jl_r.size()>=kmesh);
        for (int ik = 0; ik < kmesh; ++ik)
        {
            integrated_func[ik] = jl_r[ik] * k1_dot_k2[ik];
        }
        double temp = 0.0;

        ModuleBase::Integral::Simpson_Integral(kmesh, ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func), dk, temp);
        rs.at(ir) = temp * ModuleBase::FOUR_PI;

        const std::vector<double>& jlm1_r = jlm1.at(ir);
        const std::vector<double>& jlp1_r = jlp1.at(ir);
        assert(jlm1_r.size()>=kmesh);
        assert(jlp1_r.size()>=kmesh);
        const double fac = l / (l + 1.0);
        if (l == 0)
        {
            for (int ik = 0; ik < kmesh; ++ik)
            {
                integrated_func[ik] = jlp1_r[ik] * k1_dot_k2_dot_kpoint[ik];
            }
        }
        else
        {
            for (int ik = 0; ik < kmesh; ++ik)
            {
                integrated_func[ik] = (jlp1_r[ik] - fac * jlm1_r[ik]) * k1_dot_k2_dot_kpoint[ik];
            }
        }

        ModuleBase::Integral::Simpson_Integral(kmesh, ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func), dk, temp);
        drs.at(ir) = -ModuleBase::FOUR_PI * (l + 1) / (2.0 * l + 1) * temp;
    }

    // cal rs[0] special
    if (l > 0)
    {
        if (radials.find(0) != radials.end())
        {
            for (int ik = 0; ik < kmesh; ik++)
            {
                integrated_func[ik] = k1_dot_k2[ik] * pow(kpoint[ik], l);
            }
            double temp = 0.0;

            ModuleBase::Integral::Simpson_Integral(kmesh,
                                                   ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func),
                                                   dk,
                                                   temp);

            // PLEASE try to make dualfac function as input parameters
            // mohan note 2021-03-23
            rs.at(0) = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac(2 * l + 1) * temp;
        }
    }

    ModuleBase::timer::tick("Center2_Orb", "cal_ST_Phi12_R");

    return;
}
