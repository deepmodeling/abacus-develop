#include "../librpa_bz_sampling.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}

void test_stored_irreducible_points_keep_identity_labels()
{
    constexpr int stored_points = 8;
    require(RpaLriDetail::librpa_stored_coulomb_q_count(stored_points) == stored_points,
            "every stored Coulomb q point must be declared");
    for (int ik = 0; ik != stored_points; ++ik)
    {
        const auto index = RpaLriDetail::librpa_stored_q_index(ik, stored_points);
        require(index.coulomb_irreducible_index == ik + 1,
                "stored Coulomb q labels must be one-based and unique");
        require(index.representative_scf_index == ik + 1,
                "stored SCF representative labels must refer to their own rows");
    }
}

void test_invalid_counts_and_indices_are_rejected()
{
    bool rejected = false;
    try
    {
        (void)RpaLriDetail::librpa_stored_coulomb_q_count(0);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "an empty stored q mesh must be rejected");

    for (const int ik: {-1, 8})
    {
        rejected = false;
        try
        {
            (void)RpaLriDetail::librpa_stored_q_index(ik, 8);
        }
        catch (const std::invalid_argument&)
        {
            rejected = true;
        }
        require(rejected, "out-of-range stored q indices must be rejected");
    }
}
} // namespace

int main()
{
    test_stored_irreducible_points_keep_identity_labels();
    test_invalid_counts_and_indices_are_rejected();
    std::cout << "LibRPA BZ sampling tests passed\n";
    return 0;
}
