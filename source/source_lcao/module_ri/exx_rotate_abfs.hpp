#ifndef EXX_ROTATE_ABFS_HPP
#define EXX_ROTATE_ABFS_HPP
#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_basis/module_ao/element_basis_index-ORB.h"
#include "exx_rotate_abfs.h"

#include <cmath> // For std::erfc function used in smooth truncation
#include <vector>

template <typename Tdata>
void Moment_abfs<Tdata>::cal_multipole(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in)
{
    ModuleBase::TITLE("Rotate_abfs", "cal_multipole");
    ModuleBase::timer::start("Rotate_abfs", "cal_multipole");

    this->multipole.resize(orb_in.size());
    for (size_t T = 0; T != orb_in.size(); ++T)
    {
        this->multipole[T].resize(orb_in[T].size());
        for (size_t L = 0; L != orb_in[T].size(); ++L)
        {
            this->multipole[T][L].resize(orb_in[T][L].size());
            for (size_t N = 0; N != orb_in[T][L].size(); ++N)
            {
                const Numerical_Orbital_Lm& orb_lm = orb_in[T][L][N];
                const int nr = orb_lm.getNr();
                double* integrated_func = new double[nr];
                for (size_t ir = 0; ir != nr; ++ir)
                    integrated_func[ir] = orb_lm.getPsi(ir) * std::pow(orb_lm.getRadial(ir), 2 + L) / (2 * L + 1);

                ModuleBase::Integral::Simpson_Integral(nr, integrated_func, orb_lm.getRab(), this->multipole[T][L][N]);
            }
        }
    }

    ModuleBase::timer::end("Rotate_abfs", "cal_multipole");
}

template <typename Tdata>
void Moment_abfs<Tdata>::rotate_abfs(std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in)
{
    ModuleBase::TITLE("Rotate_abfs", "rotate_abfs");
    ModuleBase::timer::start("Rotate_abfs", "rotate_abfs");

    // construct tranformation matrix A
    for (int T = 0; T != orb_in.size(); ++T)
    {
        for (int L = 0; L != orb_in[T].size(); ++L)
        {
            for (int N = 0; N != orb_in[T][L].size(); ++N)
            {
                Numerical_Orbital_Lm& orb_lm_old = orb_in[T][L][N];
                Numerical_Orbital_Lm orb_lm_mod = orb_lm_old * 0.0;
                double square = 0.0;
                for (int N2 = 0; N2 != orb_in[T][L].size(); ++N2)
                {
                    square += this->multipole[T][L][N2] * this->multipole[T][L][N2];
                }
                double norm = std::sqrt(square);
                // change abfs
                if (N == 0)
                {
                    for (int N2 = 0; N2 != orb_in[T][L].size(); ++N2)
                        orb_lm_mod += orb_in[T][L][N2] * this->multipole[T][L][N2] / norm;
                }
                else
                {
                    for (int N2 = 0; N2 != orb_in[T][L].size(); ++N2)
                    {
                        if (N2 == N)
                        {
                            orb_lm_mod += (1.0 - this->multipole[T][L][N] * this->multipole[T][L][N2] / square)
                                          * orb_in[T][L][N2];
                        }
                        else
                        {
                            orb_lm_mod
                                += (-this->multipole[T][L][N] * this->multipole[T][L][N2] / square) * orb_in[T][L][N2];
                        }
                    }
                }
                orb_lm_old = orb_lm_mod;
                // change moment
                if (N == 0)
                {
                    this->multipole[T][L][N] = norm;
                    std::cout << "Atom type " << T << ", L " << L << ", N " << N
                              << ", multipole after rotation: " << this->multipole[T][L][N] << std::endl;
                }
                else
                {
                    this->multipole[T][L][N] = 0.0;
                }
            }
        }
    }

    ModuleBase::timer::end("Rotate_abfs", "rotate_abfs");
}

template <typename Tdata>
double Moment_abfs<Tdata>::dfact(const int& l) const
{
    double result = 1;
    for (int i = l; i > 1; i -= 2)
    {
        result *= i;
    }
    return result;
}

template <typename Tdata>
int Moment_abfs<Tdata>::factorial(const int& n) const
{
    if (n == 0)
    {
        return 1;
    }
    else if (n > 0)
    {
        return n * this->factorial(n - 1);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Moment_abfs::factorial", "n is out of range");
        return 0;
    }
}

template <typename Tdata>
double Moment_abfs<Tdata>::ln_factorial(int n) const
{
    double res = 0.0;
    for (int i = 2; i <= n; ++i)
        res += std::log(i);
    return res;
}

template <typename Tdata>
double Moment_abfs<Tdata>::cal_cl1l2(int l1, int l2) const
{
    double result = 0.0;
    // int overflow
    result = ModuleBase::FOUR_PI * std::sqrt(1.0 / ModuleBase::PI_HALF) * dfact(2 * l1 + 2 * l2 - 1) / dfact(2 * l1 - 1)
             / dfact(2 * l2 - 1);

    return result;
}

template <typename Tdata>
double Moment_abfs<Tdata>::sum_triple_Y_YLM_real(int l1,
                                                 int m1, // real m1, not index
                                                 int l2,
                                                 int m2,                         // real m2, not index
                                                 const std::vector<double>& rly, // real Y_LM(R)
                                                 const ORB_gaunt_table& MGT,
                                                 const double distance)
{
    double sum = 0.0;
    const double tiny2 = 1e-10;
    const int L = l1 + l2;
    const int idx1 = MGT.get_lm_index(l1, m1 + l1);
    const int idx2 = MGT.get_lm_index(l2, m2 + l2);

    // cyl(m1,m2) = sum_M C(l1,l2,L,m1,m2,M) Y_LM(R)
    for (int M = -L; M <= L; ++M)
    {
        const int idxL = MGT.get_lm_index(L, M + L);
        const double C = MGT.Gaunt_Coefficients(idx1, idx2, idxL);
        const double ylm_solid = rly.at(idxL);
        const double ylm_real = (distance > tiny2) ? ylm_solid / pow(distance, l1 + l2) : ylm_solid;

        if (std::abs(C) > 1e-14)
        {
            sum += C * ylm_real;
        }
    }

    return sum;
}

template <typename Tdata>
void Moment_abfs<Tdata>::cal_VR(
    const UnitCell& ucell,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in,
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>& list_r,
    const std::vector<double>& orb_cutoff,
    const double Rc,
    LRI_CV<Tdata>& cv,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut)
{
    ModuleBase::TITLE("Rotate_abfs", "cal_VR");
    ModuleBase::timer::start("Rotate_abfs", "cal_VR");
    // copy in source\module_hamilt_lcao\hamilt_lcaodft\center2_orb-orb11.cpp
    const double tiny1 = 1e-12;
    const ModuleBase::Element_Basis_Index::Range range = ModuleBase::Element_Basis_Index::construct_range(orb_in);
    // index: T: type, L: angular momentum, N: radial index, M: magnetic moment
    const ModuleBase::Element_Basis_Index::IndexLNM index = ModuleBase::Element_Basis_Index::construct_index(range);

    ORB_gaunt_table MGT0;
    int Lmax = 0;
    for (size_t T = 0; T != orb_in.size(); ++T)
    {
        Lmax = std::max(Lmax, static_cast<int>(orb_in[T].size()) - 1);
    }

    ModuleBase::TITLE("cal_Gaunt_before");
    MGT0.init_Gaunt_CH(Lmax);
    MGT0.init_Gaunt(Lmax);
    ModuleBase::TITLE("cal_Gaunt_after");

    const auto list_A0 = list_r.first;
    const auto list_A1 = list_r.second[0];

    auto& Vws = cv.Vws;
    for (size_t i0 = 0; i0 < list_A0.size(); ++i0)
    {
        const TA iat0 = list_A0[i0];
        const int T1 = ucell.iat2it[iat0];
        const size_t I1 = ucell.iat2ia[iat0];
        const auto& tauA = ucell.atoms[T1].tau[I1];
        const size_t sizeA = index[T1].count_size;
        for (size_t i1 = 0; i1 < list_A1.size(); ++i1)
        {
            const TA iat1 = list_A1[i1].first;
            const int T2 = ucell.iat2it[iat1];
            const size_t I2 = ucell.iat2ia[iat1];
            const auto& tauB = ucell.atoms[T2].tau[I2];
            const size_t sizeB = index[T2].count_size;
            const auto R = list_A1[i1].second;
            // delta_R: Angstrom
            const auto delta_R = tauB - tauA + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
            // bohr
            const double distance_true = (delta_R).norm() * ucell.lat0;
            // bohr
            const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;
            // Rcut_lcao: bohr, ucell.lat0 = 1.8897 transform angstrom to bohr
            double Rcut_lcao = orb_cutoff[T1] + orb_cutoff[T2];
            // Rcut_coul: bohr
            double Rcut_coul = std::min(cv.cal_V_Rcut(T1, T2), cv.cal_V_Rcut(T2, T2));
            if (distance < Rcut_lcao || distance >= Rcut_coul)
                continue;
            const auto JR = std::make_pair(iat1, R);
            auto tmp_tensor = RI::Tensor<Tdata>({sizeA, sizeB});
            for (int L1 = 0; L1 != orb_in[T1].size(); ++L1)
            {
                for (int L2 = 0; L2 != orb_in[T2].size(); ++L2)
                {
                    std::vector<double> rly;
                    // keep bohr
                    ModuleBase::Ylm::rl_sph_harm(L1 + L2,
                                                 (delta_R * ucell.lat0).x,
                                                 (delta_R * ucell.lat0).y,
                                                 (delta_R * ucell.lat0).z,
                                                 rly);
                    const double prefactor1 = std::pow(distance, L1 + L2 + 1);
                    for (int M1 = -L1; M1 <= L1; ++M1)
                    {
                        const int index_M1 = M1 + L1;
                        for (int M2 = -L2; M2 <= L2; ++M2)
                        {
                            const int index_M2 = M2 + L2;
                            const double prefactor = std::pow(-1, L2) * std::pow(ModuleBase::TWO_PI, 1.5) / prefactor1;
                            const double clmlm = this->cal_cl1l2(L1, L2);
                            // For real spherical harmonics, m order is: 0, 0, 1, -1, 0, 1, -1, 2, -2, ...
                            const double ylm = sum_triple_Y_YLM_real(L1, M1, L2, M2, rly, MGT0, distance);

                            // Determine N1 and N2 loop ranges based on rotate_abfs
                            // When rotate_abfs=true: only N=0 has non-zero moment, calculate only N1=0 and N2=0
                            // When rotate_abfs=false: all moments are non-zero, calculate all N1, N2
                            const int N1_max = GlobalC::exx_info.info_ri.rotate_abfs ? 1 : orb_in[T1][L1].size();
                            const int N2_max = GlobalC::exx_info.info_ri.rotate_abfs ? 1 : orb_in[T2][L2].size();

                            for (int N1 = 0; N1 != N1_max; ++N1)
                            {
                                for (int N2 = 0; N2 != N2_max; ++N2)
                                {
                                    double mom1 = this->multipole[T1][L1][N1];
                                    const double mom2 = this->multipole[T2][L2][N2];
                                    // every L has only one moment!=0 after rotation (N=0)
                                    const size_t iA = index[T1][L1][N1][index_M1];
                                    const size_t iB = index[T2][L2][N2][index_M2];

                                    const double value = prefactor * mom1 * mom2 * clmlm * ylm;

                                    // Apply smooth truncation using erfc function in log space,
                                    // exactly as implemented in FHI-aims to ensure consistency
                                    // with k-space 1-cos(qRc) truncation scheme.
                                    //
                                    // FHI-aims formula (cut_coulomb_operator.f90 lines 205, 220):
                                    //   damp = 0.5 * erfc( (ln(r) - ln(Rc)) / ln(width) )
                                    //
                                    // This gives:
                                    //   - r = Rc        : damp = 0.5
                                    //   - r = Rc/width  : damp = 0.5 * erfc(-1) ≈ 0.921
                                    //   - r = Rc*width  : damp = 0.5 * erfc(1)  ≈ 0.079
                                    //
                                    // The log-space erfc truncation and 1-cos(qRc) form a Fourier
                                    // transform pair, eliminating high-frequency oscillations
                                    // that cause aliasing in coarse k-grids.
                                    //
                                    // Reference: Spencer & Alavi, PRB 77, 193110 (2008)
                                    // Reference: FHI-aims/src/cut_coulomb_operator.f90
                                    //   - lines 71-78: analytic expression
                                    //   - lines 205, 220: implementation in log space
                                    //   - runtime_choices.f90 line 2701: width_factor = 8.0
                                    //
                                    // Default width parameter from FHI-aims (can be adjusted).
                                    // Smaller width = sharper truncation (better for coarse k-grids)
                                    // For high-lying states, stronger truncation is often needed.
                                    // width = 1.05: r = Rc*1.2 gives ~1.5% contribution
                                    // width = 1.03: r = Rc*1.1 gives ~0.5% contribution
                                    // width = 1.02: r = Rc*1.05 gives ~0.08% contribution (nearly hard)
                                    // width = 1.01: essentially hard cutoff at Rc
                                    const double width = 0.9; // Width factor - tune this for high states!

                                    double cutoff_factor = 1.0;
                                    // Uncomment to test hard cutoff (for debugging)
                                    // cutoff_factor = (distance < Rc) ? 1.0 : 0.0;
                                    // tmp_tensor(iA, iB) = value * cutoff_factor;
                                    // continue;  // Skip the erfc truncation below

                                    if (distance > 0.0 && width > 1.0)
                                    {
                                        // Log-space erfc truncation (FHI-aims implementation)
                                        cutoff_factor = 0.5 * std::erfc(std::log(distance / Rc) / std::log(width));

                                        // Debug output for high states (check if truncation is working)
                                        static int debug_count = 0;
                                        if (debug_count < 10 && distance > Rc * 0.9 && distance < Rc * 1.2)
                                        {
                                            std::cout << "DEBUG: distance=" << distance << " Rc=" << Rc
                                                      << " cutoff_factor=" << cutoff_factor << " value=" << value
                                                      << std::endl;
                                            debug_count++;
                                        }
                                    }
                                    else if (distance <= 0.0)
                                    {
                                        // At r = 0, no truncation
                                        cutoff_factor = 1.0;
                                    }
                                    else
                                    {
                                        // width <= 1.0: no smooth truncation, use hard cutoff
                                        cutoff_factor = (distance < Rc) ? 1.0 : 0.0;
                                        // const double gamma = 5.0 / Rc;
                                        // double x = gamma * distance;                           
                                        // cutoff_factor = std::erfc(x);                                       
                                        // cutoff_factor = (cutoff_factor > 0) ? cutoff_factor : 0.0;
                                    }

                                    tmp_tensor(iA, iB) = value * cutoff_factor;
                                }
                            }
                            // if (iat0 == 0 && iat0 == 0 && R[0] == 2 && R[1] == 2 && R[2] == 1 && iat0 == 0 &&
                            // iat1 == 0)
                            // {
                            //     std::cout << "iA: " << iA << ", iB: " << iB << ", L2: " << L2 << ", M2: " << M2
                            //               << ", V= " << value << ", prefactor:" << prefactor << ",mom1: " << mom1
                            //               << ", mom2: " << mom2 << ", clmlm: " << clmlm << ", ylm: " << ylm
                            //               << std::endl;
                            // }

                            // if (R[0] == 2 && R[1] == 2 && R[2] == 1 && iat0 == 0 && iat1 == 0)
                            // {
                            //     out_pure_ri_tensor("Vs_tensor.txt", tmp_tensor, 0.0);
                            // }
                            // debug
                            // std::cout << "T1: " << T1 << ", L1: " << L1 << ", T2: " << T2 << ", L2: " << L2
                            //           << ", delta_R: " << delta_R[0] << ", " << delta_R[1] << ", " << delta_R[2]
                            //           << ", distance: " << distance << ", outside V= " << value
                            //           << ", prefactor:" << prefactor << ",mom1: " << mom1 << ", mom2: " << mom2
                            //           << ", clmlm: " << clmlm << ", ylm_real: " << ylm_real << std::endl;
                        }
                    }
                }
            }
            Vws[T1][T2][delta_R] = tmp_tensor;
            // I must contain all atoms in unit cell
            auto& target_inner = Vs_cut.at(iat0);
            if (target_inner.find(JR) == target_inner.end())
            {
                target_inner.emplace(JR, tmp_tensor);
            }
            // otherwise, warning
            else
            {
                target_inner[JR] = tmp_tensor;
                const auto J = JR.first;
                const auto R = JR.second;
                // std::cout << "J=" << J << ", R= " << R[0] << ", " << R[1] << ", " << R[2] << std::endl;
                ModuleBase::WARNING_QUIT("merge_VR", "JR already exists in Vs_cut[I], which is unexpected.");
            }
        }
    }
    ModuleBase::timer::end("Rotate_abfs", "cal_VR");
}

template <typename Tdata>
void Moment_abfs<Tdata>::discard0_VR(
    const UnitCell& ucell,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in,
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>& list_r,
    const std::vector<double>& orb_cutoff,
    const double Rc,
    LRI_CV<Tdata>& cv,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut)
{
    ModuleBase::TITLE("Rotate_abfs", "discard0_VR");
    ModuleBase::timer::start("Rotate_abfs", "discard0_VR");

    const double tiny1 = 1e-12;
    const ModuleBase::Element_Basis_Index::Range range = ModuleBase::Element_Basis_Index::construct_range(orb_in);
    const ModuleBase::Element_Basis_Index::IndexLNM index = ModuleBase::Element_Basis_Index::construct_index(range);

    const auto list_A0 = list_r.first;
    const auto list_A1 = list_r.second[0];
    auto& Vws = cv.Vws;

    // Process cv.Vws - only modify tensors in the moment method range
    for (size_t i0 = 0; i0 < list_A0.size(); ++i0)
    {
        const TA iat0 = list_A0[i0];
        const int T1 = ucell.iat2it[iat0];
        const size_t I1 = ucell.iat2ia[iat0];
        const auto& tauA = ucell.atoms[T1].tau[I1];

        for (size_t i1 = 0; i1 < list_A1.size(); ++i1)
        {
            const TA iat1 = list_A1[i1].first;
            const int T2 = ucell.iat2it[iat1];
            const size_t I2 = ucell.iat2ia[iat1];
            const auto& tauB = ucell.atoms[T2].tau[I2];
            const auto R = list_A1[i1].second;

            // delta_R: Angstrom
            const auto delta_R = tauB - tauA + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
            // bohr
            const double distance_true = (delta_R).norm() * ucell.lat0;
            const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;

            double Rcut_lcao = orb_cutoff[T1] + orb_cutoff[T2];
            double Rcut_coul = std::min(cv.cal_V_Rcut(T1, T2), cv.cal_V_Rcut(T2, T2));

            // Skip if not in moment method range
            if (distance < Rcut_lcao || distance >= Rcut_coul)
                continue;

            // Check if this tensor exists in Vws
            if (Vws.find(T1) != Vws.end() && Vws[T1].find(T2) != Vws[T1].end())
            {
                auto& delta_R_map = Vws[T1][T2];
                if (delta_R_map.find(delta_R) != delta_R_map.end())
                {
                    RI::Tensor<Tdata>& tensor = delta_R_map[delta_R];

                    // Zero out elements where N1 != 0 or N2 != 0
                    for (int L1 = 0; L1 != orb_in[T1].size(); ++L1)
                    {
                        for (int L2 = 0; L2 != orb_in[T2].size(); ++L2)
                        {
                            for (int M1 = -L1; M1 <= L1; ++M1)
                            {
                                const int index_M1 = M1 + L1;
                                for (int M2 = -L2; M2 <= L2; ++M2)
                                {
                                    const int index_M2 = M2 + L2;

                                    // Set all N1 != 0 or N2 != 0 elements to zero
                                    for (int N1 = 1; N1 != orb_in[T1][L1].size(); ++N1)
                                    {
                                        const size_t iA = index[T1][L1][N1][index_M1];
                                        for (int N2 = 0; N2 != orb_in[T2][L2].size(); ++N2)
                                        {
                                            const size_t iB = index[T2][L2][N2][index_M2];
                                            tensor(iA, iB) = static_cast<Tdata>(0);
                                        }
                                    }
                                    for (int N2 = 1; N2 != orb_in[T2][L2].size(); ++N2)
                                    {
                                        const size_t iB = index[T2][L2][N2][index_M2];
                                        for (int N1 = 0; N1 != orb_in[T1][L1].size(); ++N1)
                                        {
                                            const size_t iA = index[T1][L1][N1][index_M1];
                                            tensor(iA, iB) = static_cast<Tdata>(0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Process Vs_cut - only modify tensors in the moment method range
    for (auto& iat_inner_map : Vs_cut)
    {
        const TA iat0 = iat_inner_map.first;
        const int T1 = ucell.iat2it[iat0];
        const size_t I1 = ucell.iat2ia[iat0];
        const auto& tauA = ucell.atoms[T1].tau[I1];

        for (auto& JR_tensor : iat_inner_map.second)
        {
            const TA iat1 = JR_tensor.first.first;
            const int T2 = ucell.iat2it[iat1];
            const size_t I2 = ucell.iat2ia[iat1];
            const auto& tauB = ucell.atoms[T2].tau[I2];
            const auto& R = JR_tensor.first.second;

            // delta_R: Angstrom
            const auto delta_R = tauB - tauA + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
            // bohr
            const double distance_true = (delta_R).norm() * ucell.lat0;
            const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;

            double Rcut_lcao = orb_cutoff[T1] + orb_cutoff[T2];
            double Rcut_coul = std::min(cv.cal_V_Rcut(T1, T2), cv.cal_V_Rcut(T2, T2));

            // Skip if not in moment method range
            if (distance < Rcut_lcao || distance >= Rcut_coul)
                continue;

            RI::Tensor<Tdata>& tensor = JR_tensor.second;

            // Zero out elements where N1 != 0 or N2 != 0
            for (int L1 = 0; L1 != orb_in[T1].size(); ++L1)
            {
                for (int L2 = 0; L2 != orb_in[T2].size(); ++L2)
                {
                    for (int M1 = -L1; M1 <= L1; ++M1)
                    {
                        const int index_M1 = M1 + L1;
                        for (int M2 = -L2; M2 <= L2; ++M2)
                        {
                            const int index_M2 = M2 + L2;

                            // Set all N1 != 0 or N2 != 0 elements to zero
                            for (int N1 = 1; N1 != orb_in[T1][L1].size(); ++N1)
                            {
                                const size_t iA = index[T1][L1][N1][index_M1];
                                for (int N2 = 0; N2 != orb_in[T2][L2].size(); ++N2)
                                {
                                    const size_t iB = index[T2][L2][N2][index_M2];
                                    tensor(iA, iB) = static_cast<Tdata>(0);
                                }
                            }
                            for (int N2 = 1; N2 != orb_in[T2][L2].size(); ++N2)
                            {
                                const size_t iB = index[T2][L2][N2][index_M2];
                                for (int N1 = 0; N1 != orb_in[T1][L1].size(); ++N1)
                                {
                                    const size_t iA = index[T1][L1][N1][index_M1];
                                    tensor(iA, iB) = static_cast<Tdata>(0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("Rotate_abfs", "discard0_VR");
}

template <typename Tdata>
void Moment_abfs<Tdata>::out_pure_ri_tensor(const std::string fn, RI::Tensor<double>& olp, const double threshold)
{
    std::ofstream fs;
    auto format = std::scientific;
    int prec = 15;
    fs.open(fn);
    int nr = olp.shape[0];
    int nc = olp.shape[1];
    size_t nnz = nr * nc;
    fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
    fs << "%" << std::endl;

    fs << nr << " " << nc << " " << nnz << std::endl;

    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            auto v = olp(i, j);
            if (fabs(v) > threshold)
                fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v << "\n";
        }
    }

    fs.close();
}

template <typename Tdata>
void Moment_abfs<Tdata>::out_pure_ri_tensor(const std::string fn,
                                            RI::Tensor<std::complex<double>>& olp,
                                            const double threshold)
{
    std::ofstream fs;
    auto format = std::scientific;
    int prec = 15;
    fs.open(fn);
    int nr = olp.shape[0];
    int nc = olp.shape[1];
    size_t nnz = nr * nc;
    fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
    fs << "%" << std::endl;

    fs << nr << " " << nc << " " << nnz << std::endl;

    for (int j = 0; j < nc; j++)
    {
        for (int i = 0; i < nr; i++)
        {
            auto v = olp(i, j);
            if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
                fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real()
                   << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
        }
    }

    fs.close();
}
// template <typename Tdata>
// void Moment_abfs<Tdata>::diverge_list(
//     const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>&
//     list_As_Vs, std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>&
//     list_k, std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>>& list_r,
//     const UnitCell& ucell,
//     const std::vector<double>& orb_cutoff)
// {
//     ModuleBase::TITLE("Rotate_abfs", "diverge_list");
//     ModuleBase::timer::start("Rotate_abfs", "diverge_list");

//     double Rcut_max = 0;
//     for (int T = 0; T < ucell.ntype; ++T)
//         Rcut_max = std::max(Rcut_max, orb_cutoff[T]);

//     list_k.first.clear();
//     list_k.second.clear();
//     list_r.first.clear();
//     list_r.second.clear();
//     list_k.second.resize(1);
//     list_r.second.resize(1);
//     int flag_k = 0;
//     int flag_r = 0;
//     for (size_t iA = 0; iA < list_As_Vs.first.size(); ++iA)
//     {
//         const auto& A = list_As_Vs.first[iA];
//         const size_t TA = ucell.iat2it[A];
//         const size_t IA = ucell.iat2ia[A];
//         const auto& tauA = ucell.atoms[TA].tau[IA];
//         for (const auto& BR: list_As_Vs.second[0])
//         {
//             const auto& B = BR.first;
//             const size_t TB = ucell.iat2it[B];
//             const size_t IB = ucell.iat2ia[B];
//             const auto& tauB = ucell.atoms[TB].tau[IB];
//             const auto& R = BR.second;

//             const ModuleBase::Vector3<double> tauB_shift = tauB + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
//             const ModuleBase::Vector3<double> tau_delta = (tauB_shift - tauA) * ucell.lat0;
//             const double distance = tau_delta.norm();
//             if (distance <= Rcut_max)
//             {
//                 if (std::find(list_k.first.begin(), list_k.first.end(), A) == list_k.first.end())
//                 {
//                     list_k.first.emplace_back(A);
//                 }
//                 if (std::find(list_k.second[0].begin(), list_k.second[0].end(), BR) == list_k.second[0].end())
//                 {
//                     list_k.second[0].emplace_back(BR);
//                 }
//             }
//             else
//             {
//                 if (std::find(list_r.first.begin(), list_r.first.end(), A) == list_r.first.end())
//                 {
//                     list_r.first.emplace_back(A);
//                 }
//                 if (std::find(list_r.second[0].begin(), list_r.second[0].end(), BR) == list_r.second[0].end())
//                 {
//                     list_r.second[0].emplace_back(BR);
//                 }
//             }
//         }
//     }
//     for (const auto& I: list_r.first)
//     {
//         for (const auto& JR: list_r.second[0])
//             flag_r += 1;
//     }
//     for (const auto& I: list_k.first)
//     {
//         for (const auto& JR: list_k.second[0])
//             flag_k += 1;
//     }
//     std::cout << "All No.(atom pairs)=" << flag_k + flag_r << ", " << flag_k << " inside Rcut, " << flag_r
//               << " outside Rcut" << std::endl;
//     ModuleBase::timer::end("Rotate_abfs", "diverge_list");
// }

// template <typename Tdata>
// void Moment_abfs<Tdata>::merge_list(std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut)
// {
//     ModuleBase::TITLE("Rotate_abfs", "merge_list");
//     ModuleBase::timer::start("Rotate_abfs", "merge_list");

//     for (const auto& IJRc: this->VR)
//     {
//         const TA& I = IJRc.first;
//         const auto& JRc = IJRc.second;

//         // if Vs_cut[I] does not exist, then create an empty inner map
//         auto& target_inner = Vs_cut[I];

//         for (const auto& JRc_tensor: JRc)
//         {
//             const TAC& JR = JRc_tensor.first;
//             const RI::Tensor<Tdata>& tensor = JRc_tensor.second;

//             // if JR does not exist in Vs_cut[I], then insert it
//             if (target_inner.find(JR) == target_inner.end())
//             {
//                 target_inner.emplace(JR, tensor);
//             }
//             // otherwise, warning
//             else
//             {
//                 target_inner[JR] = tensor;
//                 const auto J = JR.first;
//                 const auto R = JR.second;
//                 std::cout << "J=" << J << ", R= " << R[0] << ", " << R[1] << ", " << R[2] << std::endl;
//                 // ModuleBase::WARNING_QUIT("merge_VR", "JR already exists in Vs_cut[I], which is unexpected.");
//             }
//         }
//     }
//     ModuleBase::timer::end("Rotate_abfs", "merge_list");
// }

#endif
