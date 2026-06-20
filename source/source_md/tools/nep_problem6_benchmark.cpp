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

std::vector<std::string> parse_string_list(const std::string& text)
{
  std::vector<std::string> values;
  std::stringstream stream(text);
  std::string item;
  while (std::getline(stream, item, ',')) {
    if (!item.empty()) {
      values.push_back(item);
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

int atom_type(const int atom, const std::string& composition)
{
  if (composition == "single") {
    return 0;
  }
  if (composition == "skewed") {
    return atom % 4 == 0 ? 1 : 0;
  }
  return atom % 2;
}

SystemData make_system(const int natom, const std::string& composition)
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
        system.type.push_back(atom_type(atom, composition));
        system.position[atom] = spacing * (ix + 0.5);
        system.position[natom + atom] = spacing * (iy + 0.5);
        system.position[2 * natom + atom] = spacing * (iz + 0.5);
        ++atom;
      }
    }
  }
  return system;
}

std::vector<double> find_descriptor(NEP& nep, const SystemData& system)
{
  std::vector<double> descriptor(system.type.size() * nep.annmb.dim, 0.0);
  nep.find_descriptor(system.type, system.box, system.position, descriptor);
  return descriptor;
}

NepOutput compute(NEP& nep, const SystemData& system)
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

double time_descriptor_ms(NEP& nep, const SystemData& system, const int repeats)
{
  const auto start = std::chrono::high_resolution_clock::now();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    find_descriptor(nep, system);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count() / repeats;
}

double time_compute_ms(NEP& nep, const SystemData& system, const int repeats)
{
  const auto start = std::chrono::high_resolution_clock::now();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    compute(nep, system);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count() / repeats;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " MODEL [SIZES] [THREADS] [COMPOSITIONS] [REPEATS]\n";
    return 2;
  }

  const std::string model_file = argv[1];
  std::vector<int> sizes = parse_int_list(argc > 2 ? argv[2] : "64,216,512");
  std::vector<int> threads = parse_int_list(argc > 3 ? argv[3] : "1,2,4");
  std::vector<std::string> compositions =
    parse_string_list(argc > 4 ? argv[4] : "single,mixed,skewed");
  const int repeats = argc > 5 ? std::max(1, std::atoi(argv[5])) : 10;

  if (sizes.empty() || threads.empty() || compositions.empty()) {
    std::cerr << "Sizes, threads, and compositions must not be empty.\n";
    return 2;
  }

  std::cout
    << "atoms,composition,threads,repeats,descriptor_ms,descriptor_speedup,compute_ms,compute_speedup,"
       "max_descriptor_diff,max_potential_diff,max_force_diff,max_virial_diff\n";

  for (const int natom : sizes) {
    for (const std::string& composition : compositions) {
      const SystemData system = make_system(natom, composition);

      set_num_threads(1);
      NEP reference_nep;
      reference_nep.init_from_file(model_file, false);
      reference_nep.set_force_parallel_mode(NEP::ForceParallelMode::ThreadLocal);
      const std::vector<double> reference_descriptor = find_descriptor(reference_nep, system);
      const NepOutput reference_output = compute(reference_nep, system);
      const double serial_descriptor_ms = time_descriptor_ms(reference_nep, system, repeats);
      const double serial_compute_ms = time_compute_ms(reference_nep, system, repeats);

      for (const int num_threads : threads) {
        set_num_threads(num_threads);
        NEP nep;
        nep.init_from_file(model_file, false);
        nep.set_force_parallel_mode(NEP::ForceParallelMode::ThreadLocal);

        const std::vector<double> descriptor = find_descriptor(nep, system);
        const NepOutput output = compute(nep, system);
        const double descriptor_ms = num_threads == 1
          ? serial_descriptor_ms
          : time_descriptor_ms(nep, system, repeats);
        const double compute_ms = num_threads == 1
          ? serial_compute_ms
          : time_compute_ms(nep, system, repeats);

        std::cout << natom << ","
                  << composition << ","
                  << num_threads << ","
                  << repeats << ","
                  << descriptor_ms << ","
                  << serial_descriptor_ms / descriptor_ms << ","
                  << compute_ms << ","
                  << serial_compute_ms / compute_ms << ","
                  << max_abs_diff(reference_descriptor, descriptor) << ","
                  << max_abs_diff(reference_output.potential, output.potential) << ","
                  << max_abs_diff(reference_output.force, output.force) << ","
                  << max_abs_diff(reference_output.virial, output.virial) << "\n";
      }
    }
  }

  return 0;
}
