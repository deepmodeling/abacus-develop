#include "source_hamilt/module_vdw/vdwd3_data.h"
#include "source_hamilt/module_vdw/vdwd3_evaluator.h"
#include "source_hamilt/module_vdw/vdwd3_parameters.h"
#include "source_hamilt/module_vdw/vdw_xcname.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace
{

using vdw::d3::Cutoffs;
using vdw::d3::Damping;
using vdw::d3::Parameters;
using vdw::d3::Result;
using vdw::d3::Structure;
using vdw::d3::Vec3;

double& coordinate(Vec3& vector, int component)
{
    if (component == 0)
    {
        return vector.x;
    }
    if (component == 1)
    {
        return vector.y;
    }
    return vector.z;
}

double coordinate(const Vec3& vector, int component)
{
    if (component == 0)
    {
        return vector.x;
    }
    if (component == 1)
    {
        return vector.y;
    }
    return vector.z;
}

double energy(const Structure& structure, const Parameters& parameters, const Cutoffs& cutoffs)
{
    Result result;
    std::string error;
    EXPECT_TRUE(vdw::d3::evaluate(structure, parameters, cutoffs, false, result, error)) << error;
    return result.energy;
}

Structure mindless01()
{
    Structure structure;
    structure.atomic_numbers = {11, 1, 8, 1, 9, 1, 1, 8, 7, 1, 1, 17, 5, 5, 7, 13};
    structure.positions = {
        {-1.85528263484662, 3.58670515364616, -2.41763729306344},
        {4.40178023537845, 0.02338844412653, -4.95457749372945},
        {-2.98706033463438, 4.76252065456814, 1.27043301573532},
        {0.79980886075526, 1.41103455609189, -5.04655321620119},
        {-4.20647469409936, 1.84275767548460, 4.55038084858449},
        {-3.54356121843970, -3.18835665176557, 1.46240021785588},
        {2.70032160109941, 1.06818452504054, -1.73234650374438},
        {3.73114088824361, -2.07001543363453, 2.23160937604731},
        {-1.75306819230397, 0.35951417150421, 1.05323406177129},
        {5.41755788583825, -1.57881830078929, 1.75394002750038},
        {-2.23462868255966, -2.13856505054269, 4.10922285746451},
        {1.01565866207568, -3.21952154552768, -3.36050963020778},
        {2.42119255723593, 0.26626435093114, -3.91862474360560},
        {-3.02526098819107, 2.53667889095925, 2.31664984740423},
        {-2.00438948664892, -2.29235136977220, 2.19782807357059},
        {1.12226554109716, -1.36942007032045, 0.48455055461782},
    };
    return structure;
}

Structure mindless17()
{
    Structure structure;
    structure.atomic_numbers = {12, 5, 1, 6, 14, 5, 9, 9, 1, 1, 9, 1, 1, 8, 1, 1};
    structure.positions = {
        {2.68460861953273, -1.46442870252660, 1.67856153340039},
        {-0.18053173984261, -3.96944005202894, -0.80661420861966},
        {-0.58984435587144, -5.98708837691298, -1.71743303841248},
        {1.07719032312989, -1.76303169348468, -2.22950669619910},
        {0.24421210314383, 1.76887707010948, -4.26437232785349},
        {-1.64589939447879, -0.90586156069035, -2.45065764617243},
        {0.71078922997216, 0.88479326169415, 3.64517056373104},
        {5.77984092031648, -2.28654280723977, 2.67325530102119},
        {-3.50012299703202, -1.32272929604435, -1.24933974739626},
        {-0.52632555428452, -3.77091486672575, 1.47071730349011},
        {-0.27536787387566, 4.17239735599888, -2.08796280751818},
        {-1.56874517509641, 4.33699624911296, 0.52992881707066},
        {-1.79332867160074, 5.69507449430072, 3.22335048552804},
        {-2.07387284201988, 4.10380585322630, 2.35332246894739},
        {-0.64011933245457, 2.42620166307914, 3.12632079492992},
        {2.29751674046149, -1.91810859186826, -3.89474079594712},
    };
    return structure;
}

Structure mindless09()
{
    Structure structure;
    structure.atomic_numbers = {1, 1, 1, 1, 3, 1, 6, 5, 1, 1, 14, 1, 17, 9, 1, 5};
    structure.positions = {
        {3.97360649552839, 1.71723751297383, -0.51862929250676},
        {0.16903666216522, 1.73154352333176, -0.40099024352959},
        {-3.94463844105182, -1.24346369608005, 0.09565841726334},
        {2.21647168119803, 4.10625979391554, 2.61391340002321},
        {-0.04488993380842, -2.16288302687041, 4.48488595610432},
        {3.52287141817194, -0.90500888687059, -5.00916337263077},
        {1.95336082370762, -0.83849036872324, -3.65515970516029},
        {2.05706981818495, 1.70095588601056, -2.06303335904159},
        {-6.40097100472159, -1.71072935987273, 3.14621771036234},
        {2.04751538182937, -2.55691868000982, -2.49926722310562},
        {2.03251078714394, 1.35094356516468, 2.02150308748654},
        {0.20477572129201, -0.93291693232462, -4.76431390827476},
        {-2.67673272939098, 1.40764602033672, 4.10347165469140},
        {-2.75901984658887, -3.73954809548334, 3.19373273207227},
        {1.96938102642596, 3.74070925169244, -3.03185101883736},
        {-4.32034786008576, -1.66533650719069, 2.28302516508337},
    };
    return structure;
}

Structure actinides()
{
    Structure structure;
    structure.atomic_numbers = {87, 88, 89, 90, 91, 92, 93, 94, 95,
                                96, 97, 98, 99, 100, 101, 102, 103};
    structure.positions = {
        {0.98692316414074, 6.12727238368797, -6.67861597188102},
        {3.63898862390869, 5.12109301182962, 3.01908613326278},
        {5.14503571563551, -3.97172984617710, 3.82011791828867},
        {6.71986847575494, 1.71382138402812, 3.92749159076307},
        {4.13783589704826, -2.10695793491818, 0.19753203068899},
        {8.97685097698326, -3.08813636191844, -4.45568615593938},
        {12.5486412940776, -1.77128765259458, 0.59261498922861},
        {7.82051475868325, -3.97159756604558, -0.53637703616916},
        {-0.43444574624893, -1.69696511583960, -1.65898182093050},
        {-4.71270645149099, -0.11534827468942, 2.84863373521297},
        {-2.52061680335614, 1.82937752749537, -2.10366982879172},
        {0.13551154616576, 7.99805359235043, -1.55508522619903},
        {3.91594542499717, -1.72975169129597, -5.07944366756113},
        {-1.03393930231679, 4.69307230054046, 0.02656940927472},
        {6.20675384557240, 4.24490721493632, -0.71004195169885},
        {7.04586341131562, 5.20053667939076, -7.51972863675876},
        {2.01082807362334, 1.34838807211157, -4.70482633508447},
    };
    return structure;
}

Structure x23_acetic()
{
    Structure structure;
    structure.atomic_numbers = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        6, 6, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8,
    };
    structure.positions = {
        {7.90377778530531, 7.12936818903081, 9.56314593193451},
        {16.95915441875529, 0.25435708066996, 4.16571135682655},
        {20.33533837362116, 3.94612522923480, 9.56314593193451},
        {4.52759383043945, 3.43741106789489, 4.16571135682655},
        {9.30463145470085, 3.53964522884767, -5.36209670431652},
        {15.55848972193083, -3.53964522884767, 10.83020702100736},
        {21.73600307044562, 7.53584818941794, 5.43277244589940},
        {11.68606379541616, 3.34481450806704, 3.04226942177239},
        {13.17705738121553, -3.34481450806704, 8.43970399688035},
        {24.11762438373200, 15.11421520732826, 3.04226942177239},
        {0.74549679289969, -0.34714261306888, 8.43970399688035},
        {11.18925490605411, 6.24818909009767, 4.59108861432140},
        {13.67386627057758, 1.13553617960310, 9.98833421685828},
        {23.62081549436995, 4.82730432816794, 4.59108861432140},
        {1.24230568226174, 9.93995723866251, -0.80634596078656},
        {3.12692913361499, -0.15212291971717, 10.83020702100736},
        {8.32726531708940, 5.25608309194219, 1.90276481817667},
        {16.53585585954229, 2.12764217775858, 7.30019939328463},
        {20.75863693283416, 5.81941032632342, 1.90276481817667},
        {4.10429527122645, 1.56412597080627, 7.30019939328463},
        {10.22719554669991, 4.55253820982165, 3.88773270477194},
        {14.63592562993178, 2.83099808730804, 9.28516727987990},
        {22.65875613501575, 6.52276623587289, 3.88773270477194},
        {2.20436504161594, 0.86077006125680, 9.28516727987990},
        {9.27212817247557, 6.62140991797520, 0.03533787079144},
        {15.59099300415612, 0.76212637915449, 5.43277244589940},
        {21.70368876079141, 4.45389452771933, 0.03533787079144},
        {3.15943241584028, 2.92964176941036, 5.43277244589940},
        {6.07055487328507, 4.63171771710301, 1.93885857925242},
        {18.79237733077554, 2.75181858002668, 7.33629315436038},
        {18.50211546160091, 6.44377570116260, 1.93885857925242},
        {6.36100571503078, 0.93994956853816, 7.33629315436038},
    };
    structure.lattice = {{{24.86304558760325, 0.0, 0.0},
                          {0.0, 7.38360999643241, 0.0},
                          {0.0, 0.0, 10.79482606446971}}};
    structure.periodic = {{true, true, true}};
    return structure;
}

} // namespace

TEST(VdwXcName, NormalizesAbacusAndLibxcSpellings)
{
    EXPECT_EQ(vdw::normalize_xc_name("PBE"), "pbe");
    EXPECT_EQ(vdw::normalize_xc_name(" XC_HYB_GGA_XC_PBEH "), "hyb_gga_xc_pbeh");
    EXPECT_EQ(vdw::normalize_xc_name("GGA_X_PBE+GGA_C_PBE"),
              "gga_x_pbe:gga_c_pbe");
    EXPECT_EQ(vdw::normalize_xc_name("XC_GGA_X_PBE + XC_GGA_C_PBE"),
              "gga_x_pbe:gga_c_pbe");
    EXPECT_EQ(vdw::normalize_xc_name("XC_GGA_X_PBE_R:XC_GGA_C_PBE"),
              "gga_x_pbe_r:gga_c_pbe");
}

TEST(D3Parameters, PbeAndLibxcAliases)
{
    Parameters zero;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Zero, zero, canonical));
    EXPECT_EQ(canonical, "pbe");
    EXPECT_DOUBLE_EQ(zero.s6, 1.0);
    EXPECT_DOUBLE_EQ(zero.s8, 0.722);
    EXPECT_DOUBLE_EQ(zero.rs6, 1.217);
    EXPECT_DOUBLE_EQ(zero.rs8, 1.0);

    Parameters rational;
    ASSERT_TRUE(vdw::d3::lookup_parameters("GGA_X_PBE+GGA_C_PBE",
                                           Damping::Rational,
                                           rational,
                                           canonical));
    EXPECT_EQ(canonical, "pbe");
    EXPECT_DOUBLE_EQ(rational.s8, 0.7875);
    EXPECT_DOUBLE_EQ(rational.a1, 0.4289);
    EXPECT_DOUBLE_EQ(rational.a2, 4.4407);

    ASSERT_TRUE(vdw::d3::lookup_parameters("HSE", Damping::Rational, rational, canonical));
    EXPECT_EQ(canonical, "hse06");
    EXPECT_DOUBLE_EQ(rational.a1, 0.383);

    ASSERT_TRUE(vdw::d3::lookup_parameters("XC_GGA_X_B88+XC_GGA_C_OP_B88",
                                           Damping::Rational,
                                           rational,
                                           canonical));
    EXPECT_EQ(canonical, "bop");
    EXPECT_DOUBLE_EQ(rational.a1, 0.487);

    ASSERT_TRUE(vdw::d3::lookup_parameters("XC_GGA_X_PBE_R+XC_GGA_C_PBE",
                                           Damping::Rational,
                                           rational,
                                           canonical));
    EXPECT_EQ(canonical, "revpbe");

    // Legacy ABACUS LibXC spellings: the numerical parameter set is stored
    // under `wb97x` by s-dftd3, but only for the corresponding damping form.
    ASSERT_TRUE(vdw::d3::lookup_parameters("XC_HYB_GGA_XC_WB97X_V",
                                           Damping::Rational,
                                           rational,
                                           canonical));
    EXPECT_EQ(canonical, "wb97xv");
    EXPECT_DOUBLE_EQ(rational.s8, 0.2641);
    EXPECT_DOUBLE_EQ(rational.a1, 0.0);
    EXPECT_DOUBLE_EQ(rational.a2, 5.4959);

    ASSERT_TRUE(vdw::d3::lookup_parameters("XC_HYB_GGA_XC_WB97X_D3",
                                           Damping::Zero,
                                           zero,
                                           canonical));
    EXPECT_EQ(canonical, "wb97xd3");
    EXPECT_DOUBLE_EQ(zero.s8, 1.0);
    EXPECT_DOUBLE_EQ(zero.rs6, 1.281);
    EXPECT_DOUBLE_EQ(zero.rs8, 1.094);

    EXPECT_FALSE(vdw::d3::lookup_parameters("XC_HYB_GGA_XC_WB97X_V",
                                            Damping::Zero,
                                            zero,
                                            canonical));
    EXPECT_FALSE(vdw::d3::lookup_parameters("XC_HYB_GGA_XC_WB97X_D3",
                                            Damping::Rational,
                                            rational,
                                            canonical));

    EXPECT_FALSE(vdw::d3::lookup_parameters("definitely-not-an-xc",
                                            Damping::Rational,
                                            rational,
                                            canonical));
}

TEST(D3Data, ReferencePackingAndSymmetry)
{
    EXPECT_EQ(vdw::d3::data::max_element, 103);
    EXPECT_EQ(vdw::d3::data::max_reference, 7);

    int maximum_reference_count = 0;
    for (int atomic_number_i = 1; atomic_number_i <= vdw::d3::data::max_element; ++atomic_number_i)
    {
        const int references_i = vdw::d3::data::reference_count(atomic_number_i);
        ASSERT_GE(references_i, 1);
        ASSERT_LE(references_i, vdw::d3::data::max_reference);
        maximum_reference_count = std::max(maximum_reference_count, references_i);
        EXPECT_TRUE(std::isfinite(vdw::d3::data::covalent_radius(atomic_number_i)));
        EXPECT_TRUE(std::isfinite(vdw::d3::data::r4r2(atomic_number_i)));
        for (int reference_i = 0; reference_i < references_i; ++reference_i)
        {
            EXPECT_TRUE(std::isfinite(vdw::d3::data::reference_cn(atomic_number_i, reference_i)));
        }

        for (int atomic_number_j = 1; atomic_number_j <= atomic_number_i; ++atomic_number_j)
        {
            const int references_j = vdw::d3::data::reference_count(atomic_number_j);
            EXPECT_DOUBLE_EQ(vdw::d3::data::vdw_radius(atomic_number_i, atomic_number_j),
                             vdw::d3::data::vdw_radius(atomic_number_j, atomic_number_i));
            for (int reference_i = 0; reference_i < references_i; ++reference_i)
            {
                for (int reference_j = 0; reference_j < references_j; ++reference_j)
                {
                    const double c6 = vdw::d3::data::reference_c6(atomic_number_i,
                                                                 reference_i,
                                                                 atomic_number_j,
                                                                 reference_j);
                    EXPECT_TRUE(std::isfinite(c6));
                    EXPECT_GE(c6, 0.0);
                    EXPECT_DOUBLE_EQ(c6,
                                     vdw::d3::data::reference_c6(atomic_number_j,
                                                                 reference_j,
                                                                 atomic_number_i,
                                                                 reference_i));
                }
            }
        }
    }
    EXPECT_EQ(maximum_reference_count, vdw::d3::data::max_reference);
    EXPECT_EQ(vdw::d3::data::reference_count(103), 7);
}

TEST(D3Evaluator, MatchesSdftd3Mindless01)
{
    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));
    parameters.s9 = 0.0;

    // Exact value from s-dftd3 v1.5.0 test/unit/test_dftd3.f90.
    EXPECT_DOUBLE_EQ(energy(mindless01(), parameters, Cutoffs()),
                     -1.7882220155186028e-2);
}

TEST(D3Evaluator, MatchesSdftd3Actinides)
{
    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));
    parameters.s9 = 0.0;

    // Exact value from s-dftd3 v1.5.0 test/unit/test_dftd3.f90.
    EXPECT_DOUBLE_EQ(energy(actinides(), parameters, Cutoffs()),
                     -1.4131143363689097e-1);
}

TEST(D3Evaluator, MatchesSdftd3PeriodicX23Acetic)
{
    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));
    parameters.s9 = 0.0;

    Cutoffs cutoffs;
    cutoffs.cn = 30.0;
    cutoffs.disp2 = 60.0;
    cutoffs.disp3 = 15.0;

    // Exact value from s-dftd3 v1.5.0 test/unit/test_periodic_3d.f90.
    EXPECT_DOUBLE_EQ(energy(x23_acetic(), parameters, cutoffs),
                     -6.6732836815486210e-2);
}

TEST(D3Evaluator, MatchesSdftd3AtmMindless17)
{
    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("dsd-blyp",
                                           Damping::Rational,
                                           parameters,
                                           canonical));
    parameters.s9 = 1.0;

    // Exact value from s-dftd3 v1.5.0 test/unit/test_dftd3.f90.
    EXPECT_DOUBLE_EQ(energy(mindless17(), parameters, Cutoffs()),
                     -1.3592755832923201e-2);
}

TEST(D3Evaluator, MatchesSdftd3ZeroDampingMindless09)
{
    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("rpbe", Damping::Zero, parameters, canonical));
    parameters.s9 = 0.0;

    // Exact value from s-dftd3 v1.5.0 test/unit/test_dftd3.f90.
    EXPECT_DOUBLE_EQ(energy(mindless09(), parameters, Cutoffs()),
                     -2.0178760785797962e-2);
}

TEST(D3Evaluator, SmoothCutoffGradientAndVirial)
{
    Structure structure;
    structure.atomic_numbers = {6, 8, 7};
    structure.positions = {{0.0, 0.0, 0.0}, {5.5, 0.4, 0.2}, {2.7, 4.8, -0.3}};

    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));
    parameters.s9 = 1.0;

    Cutoffs cutoffs;
    cutoffs.disp2 = 7.0;
    cutoffs.disp3 = 7.0;
    cutoffs.cn = 8.0;
    cutoffs.width2 = 2.0;
    cutoffs.width3 = 2.0;

    Result analytic;
    std::string error;
    ASSERT_TRUE(vdw::d3::evaluate(structure, parameters, cutoffs, true, analytic, error)) << error;

    const double step = 1.0e-5;
    double maximum_gradient_error = 0.0;
    for (std::size_t atom = 0; atom < structure.positions.size(); ++atom)
    {
        for (int component = 0; component < 3; ++component)
        {
            Structure plus = structure;
            Structure minus = structure;
            coordinate(plus.positions[atom], component) += step;
            coordinate(minus.positions[atom], component) -= step;
            const double numerical = (energy(plus, parameters, cutoffs)
                                      - energy(minus, parameters, cutoffs))
                                     / (2.0 * step);
            maximum_gradient_error = std::max(
                maximum_gradient_error,
                std::abs(numerical - coordinate(analytic.gradient[atom], component)));
        }
    }
    EXPECT_LT(maximum_gradient_error, 1.0e-9);

    double maximum_virial_error = 0.0;
    for (int strain_row = 0; strain_row < 3; ++strain_row)
    {
        for (int strain_column = 0; strain_column < 3; ++strain_column)
        {
            Structure plus = structure;
            Structure minus = structure;
            for (std::size_t atom = 0; atom < structure.positions.size(); ++atom)
            {
                const double displacement = step * coordinate(structure.positions[atom], strain_column);
                coordinate(plus.positions[atom], strain_row) += displacement;
                coordinate(minus.positions[atom], strain_row) -= displacement;
            }
            const double numerical = (energy(plus, parameters, cutoffs)
                                      - energy(minus, parameters, cutoffs))
                                     / (2.0 * step);
            maximum_virial_error = std::max(
                maximum_virial_error,
                std::abs(numerical - analytic.virial.value[strain_row][strain_column]));
        }
    }
    EXPECT_LT(maximum_virial_error, 1.0e-8);
}

TEST(D3Evaluator, PeriodicSmoothCutoffVirial)
{
    Structure structure;
    structure.atomic_numbers = {6, 8, 7};
    structure.positions = {{1.1, 1.4, 1.7}, {5.0, 1.8, 2.1}, {2.6, 5.4, 4.9}};
    structure.lattice = {{{8.6, 0.2, 0.1}, {0.3, 8.2, 0.4}, {0.2, 0.1, 8.9}}};
    structure.periodic = {{true, true, true}};

    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));
    parameters.s9 = 1.0;

    Cutoffs cutoffs;
    cutoffs.disp2 = 7.0;
    cutoffs.disp3 = 7.0;
    cutoffs.cn = 7.5;
    cutoffs.width2 = 1.5;
    cutoffs.width3 = 1.5;

    Result analytic;
    std::string error;
    ASSERT_TRUE(vdw::d3::evaluate(structure, parameters, cutoffs, true, analytic, error)) << error;

    const double step = 1.0e-5;
    double maximum_gradient_error = 0.0;
    for (std::size_t atom = 0; atom < structure.positions.size(); ++atom)
    {
        for (int component = 0; component < 3; ++component)
        {
            Structure plus = structure;
            Structure minus = structure;
            coordinate(plus.positions[atom], component) += step;
            coordinate(minus.positions[atom], component) -= step;
            const double numerical = (energy(plus, parameters, cutoffs)
                                      - energy(minus, parameters, cutoffs))
                                     / (2.0 * step);
            maximum_gradient_error = std::max(
                maximum_gradient_error,
                std::abs(numerical - coordinate(analytic.gradient[atom], component)));
        }
    }
    EXPECT_LT(maximum_gradient_error, 1.0e-8);

    double maximum_virial_error = 0.0;
    for (int strain_row = 0; strain_row < 3; ++strain_row)
    {
        for (int strain_column = 0; strain_column < 3; ++strain_column)
        {
            Structure plus = structure;
            Structure minus = structure;
            for (std::size_t atom = 0; atom < structure.positions.size(); ++atom)
            {
                const double displacement = step * coordinate(structure.positions[atom], strain_column);
                coordinate(plus.positions[atom], strain_row) += displacement;
                coordinate(minus.positions[atom], strain_row) -= displacement;
            }
            for (int vector = 0; vector < 3; ++vector)
            {
                const double displacement = step * coordinate(structure.lattice[vector], strain_column);
                coordinate(plus.lattice[vector], strain_row) += displacement;
                coordinate(minus.lattice[vector], strain_row) -= displacement;
            }
            const double numerical = (energy(plus, parameters, cutoffs)
                                      - energy(minus, parameters, cutoffs))
                                     / (2.0 * step);
            maximum_virial_error = std::max(
                maximum_virial_error,
                std::abs(numerical - analytic.virial.value[strain_row][strain_column]));
        }
    }
    EXPECT_LT(maximum_virial_error, 1.0e-8);
}

TEST(D3Evaluator, SharpCutoffCompatibilityAndWidthValidation)
{
    Structure structure;
    structure.atomic_numbers = {6, 8};
    structure.positions = {{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}};

    Parameters parameters;
    std::string canonical;
    ASSERT_TRUE(vdw::d3::lookup_parameters("pbe", Damping::Rational, parameters, canonical));

    Cutoffs sharp;
    sharp.disp2 = 7.0;
    sharp.disp3 = 7.0;
    sharp.cn = 8.0;
    sharp.width2 = 0.0;
    sharp.width3 = 0.0;

    // s-dftd3 treats width == cutoff as the legacy sharp-cutoff branch.
    Cutoffs equal_width = sharp;
    equal_width.width2 = equal_width.disp2;
    equal_width.width3 = equal_width.disp3;
    EXPECT_DOUBLE_EQ(energy(structure, parameters, sharp),
                     energy(structure, parameters, equal_width));

    Cutoffs invalid = sharp;
    invalid.width2 = invalid.disp2 + 1.0;
    Result result;
    std::string error;
    EXPECT_FALSE(vdw::d3::evaluate(structure, parameters, invalid, false, result, error));
    EXPECT_NE(error.find("0 <= width <= cutoff"), std::string::npos);
}
