#include "../potential/ml/nep/nep_cpu.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{

struct SystemData {
  std::vector<int> type;
  std::vector<double> box;
  std::vector<double> position;
};

struct NepOutput {
  std::vector<double> potential;
  std::vector<double> force;
  std::vector<double> virial;
};

std::vector<int> parse_int_list(const std::string& text)
{
  std::vector<int> values;
  std::stringstream stream(text);
  std::string item;

  while (std::getline(stream, item, ',')) {
    const int value = std::atoi(item.c_str());
    if (value > 0) {
      values.push_back(value);
    }
  }

  return values;
}

void set_num_threads(const int threads)
{
#if defined(_OPENMP)
  omp_set_num_threads(threads);
#else
  (void)threads;
#endif
}

SystemData make_system(const int natom)
{
  const int ngrid = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(natom))));
  const double spacing = 2.7;
  const double box_length = spacing * ngrid;

  SystemData system;
  system.type.reserve(natom);
  system.box = {
    box_length, 0.0, 0.0,
    0.0, box_length, 0.0,
    0.0, 0.0, box_length};
  system.position.assign(3 * natom, 0.0);

  int atom = 0;
  for (int ix = 0; ix < ngrid && atom < natom; ++ix) {
    for (int iy = 0; iy < ngrid && atom < natom; ++iy) {
      for (int iz = 0; iz < ngrid && atom < natom; ++iz) {
        system.type.push_back(atom % 2);
        system.position[atom] = spacing * (ix + 0.5);
        system.position[natom + atom] = spacing * (iy + 0.5);
        system.position[2 * natom + atom] = spacing * (iz + 0.5);
        ++atom;
      }
    }
  }

  return system;
}

NepOutput compute_once(NEP& nep, const SystemData& system)
{
  NepOutput output;
  output.potential.assign(system.type.size(), 0.0);
  output.force.assign(system.type.size() * 3, 0.0);
  output.virial.assign(system.type.size() * 9, 0.0);

  nep.compute(system.type, system.box, system.position, output.potential, output.force, output.virial);
  return output;
}

double max_abs_diff(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
  double result = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    result = std::max(result, std::abs(lhs[i] - rhs[i]));
  }
  return result;
}

double time_compute_ms(NEP& nep, const SystemData& system, const int repeats)
{
  const auto start = std::chrono::high_resolution_clock::now();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    compute_once(nep, system);
  }
  const auto end = std::chrono::high_resolution_clock::now();

  const std::chrono::duration<double, std::milli> elapsed = end - start;
  return elapsed.count() / repeats;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " MODEL_FILE [SIZES=64,216,512] [THREADS=1,2,4] [REPEATS=10]\n";
    return 2;
  }

  const std::string model_file = argv[1];
  std::vector<int> sizes = parse_int_list(argc > 2 ? argv[2] : "64,216,512");
  std::vector<int> threads = parse_int_list(argc > 3 ? argv[3] : "1,2,4");
  const int repeats = argc > 4 ? std::max(1, std::atoi(argv[4])) : 10;

  if (sizes.empty()) {
    sizes.push_back(64);
  }
  if (threads.empty()) {
    threads.push_back(1);
  }

  std::cout << "atoms,threads,repeats,avg_ms,speedup,max_potential_diff,max_force_diff,max_virial_diff\n";

  for (const int natom : sizes) {
    const SystemData system = make_system(natom);

    set_num_threads(1);
    NEP reference_nep;
    reference_nep.init_from_file(model_file, false);
    const NepOutput reference = compute_once(reference_nep, system);

    double serial_ms = 0.0;
    std::vector<double> avg_times;

    for (const int num_threads : threads) {
      set_num_threads(num_threads);
      NEP nep;
      nep.init_from_file(model_file, false);

      const NepOutput output = compute_once(nep, system);
      const double avg_ms = time_compute_ms(nep, system, repeats);

      if (num_threads == 1 || serial_ms == 0.0) {
        serial_ms = avg_ms;
      }

      const double speedup = serial_ms / avg_ms;
      std::cout << natom << ","
                << num_threads << ","
                << repeats << ","
                << avg_ms << ","
                << speedup << ","
                << max_abs_diff(reference.potential, output.potential) << ","
                << max_abs_diff(reference.force, output.force) << ","
                << max_abs_diff(reference.virial, output.virial) << "\n";
    }
  }

  return 0;
}
