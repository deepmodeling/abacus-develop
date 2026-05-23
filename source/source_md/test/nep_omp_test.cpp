#include "gtest/gtest.h"

#include "../potential/ml/nep/nep_cpu.h"

#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{

struct NepOutput {
  std::vector<double> potential;
  std::vector<double> force;
  std::vector<double> virial;
};

void set_omp_threads(const int num_threads)
{
#if defined(_OPENMP)
  omp_set_num_threads(num_threads);
#else
  (void)num_threads;
#endif
}

NepOutput compute_nep(const int num_threads)
{
  set_omp_threads(num_threads);

  NEP nep;
  nep.init_from_file(NEP_TEST_MODEL_FILE, false);

  const std::vector<int> type = {0, 1, 0, 1, 0, 1};
  const std::vector<double> box = {
    10.0, 0.0, 0.0,
    0.0, 10.0, 0.0,
    0.0, 0.0, 10.0};

  const std::vector<double> position = {
    1.0, 3.1, 5.2, 6.9, 2.4, 7.3,
    1.1, 2.8, 5.0, 7.1, 6.5, 3.7,
    0.8, 3.5, 4.9, 6.8, 7.2, 2.6};

  NepOutput output;
  output.potential.assign(type.size(), 0.0);
  output.force.assign(type.size() * 3, 0.0);
  output.virial.assign(type.size() * 9, 0.0);

  nep.compute(type, box, position, output.potential, output.force, output.virial);
  return output;
}

void expect_vectors_near(
  const std::vector<double>& reference,
  const std::vector<double>& actual,
  const double tolerance)
{
  ASSERT_EQ(reference.size(), actual.size());

  for (std::size_t i = 0; i < reference.size(); ++i) {
    EXPECT_NEAR(reference[i], actual[i], tolerance) << "index = " << i;
  }
}

} // namespace

TEST(NEPOpenMPTest, ForcePathMatchesSingleThread)
{
  const NepOutput serial = compute_nep(1);
  const NepOutput parallel = compute_nep(4);

  expect_vectors_near(serial.potential, parallel.potential, 1.0e-12);
  expect_vectors_near(serial.force, parallel.force, 1.0e-10);
  expect_vectors_near(serial.virial, parallel.virial, 1.0e-10);
}
