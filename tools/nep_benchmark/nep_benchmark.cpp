#include "nep.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

struct Result {
  double energy = 0.0;
  std::vector<double> force;
  std::vector<double> virial;
};

Result run_nep(NEP& nep, const std::vector<int>& type, const std::vector<double>& box,
               const std::vector<double>& coord)
{
  const std::size_t N = type.size();
  Result out;
  out.force.assign(N * 3, 0.0);
  out.virial.assign(N * 9, 0.0);
  std::vector<double> potential(N, 0.0);
  nep.compute(type, box, coord, potential, out.force, out.virial);
  out.energy = 0.0;
  for (double e : potential) {
    out.energy += e;
  }
  return out;
}

double max_force_diff(const Result& a, const Result& b)
{
  double max_diff = 0.0;
  for (std::size_t i = 0; i < a.force.size(); ++i) {
    max_diff = std::max(max_diff, std::abs(a.force[i] - b.force[i]));
  }
  return max_diff;
}

void build_test_system(const NEP& nep, const int natom, std::vector<int>& type,
                       std::vector<double>& box, std::vector<double>& coord)
{
  const int ntype = static_cast<int>(nep.paramb.num_types);
  type.resize(natom);
  coord.assign(natom * 3, 0.0);

  const int nx = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(natom))));
  const int ny = nx;
  const int nz = nx;
  const double a = 5.2;
  box = {a * nx, 0.0, 0.0, 0.0, a * ny, 0.0, 0.0, 0.0, a * nz};

  int iat = 0;
  for (int ix = 0; ix < nx && iat < natom; ++ix) {
    for (int iy = 0; iy < ny && iat < natom; ++iy) {
      for (int iz = 0; iz < nz && iat < natom; ++iz) {
        type[iat] = (ntype > 1) ? (iat % ntype) : 0;
        coord[iat] = (ix + 0.1 * (type[iat] + 1)) * a;
        coord[iat + natom] = (iy + 0.2 * (type[iat] + 1)) * a;
        coord[iat + 2 * natom] = (iz + 0.3 * (type[iat] + 1)) * a;
        ++iat;
      }
    }
  }
}

} // namespace

int main(int argc, char** argv)
{
  std::string model = "tests/PP_ORB/nep_hfo2.txt";
  int natom = 64;
  int repeat = 10;
  bool verify = false;
  bool perf = false;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--model" && i + 1 < argc) {
      model = argv[++i];
    } else if (arg == "--natom" && i + 1 < argc) {
      natom = std::atoi(argv[++i]);
    } else if (arg == "--repeat" && i + 1 < argc) {
      repeat = std::atoi(argv[++i]);
    } else if (arg == "--verify") {
      verify = true;
    } else if (arg == "--perf") {
      perf = true;
    }
  }

  if (!verify && !perf) {
    verify = true;
    perf = true;
  }

  NEP nep(model);
  std::vector<int> type;
  std::vector<double> box;
  std::vector<double> coord;
  build_test_system(nep, natom, type, box, coord);

  if (verify) {
#ifdef _OPENMP
    const int prev_threads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    const Result ref = run_nep(nep, type, box, coord);
#ifdef _OPENMP
    omp_set_num_threads(prev_threads);
#endif
    const Result cur = run_nep(nep, type, box, coord);
    const double e_diff = std::abs(ref.energy - cur.energy);
    const double f_diff = max_force_diff(ref, cur);
    std::cout << "energy=" << cur.energy << " dE=" << e_diff << " max_dF=" << f_diff << std::endl;
    if (e_diff > 1e-8 || f_diff > 1e-8) {
      std::cerr << "Correctness check failed: dE=" << e_diff << " max_dF=" << f_diff << std::endl;
      return 1;
    }
  }

  if (perf) {
    run_nep(nep, type, box, coord);
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < repeat; ++i) {
      run_nep(nep, type, box, coord);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double sec =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    std::cout << "natom=" << natom << " repeat=" << repeat << " total_s=" << sec
              << " per_call_ms=" << (sec * 1000.0 / repeat) << std::endl;
  }

  return 0;
}
