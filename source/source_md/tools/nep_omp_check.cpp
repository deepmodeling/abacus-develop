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

struct NepOutput {
  std::vector<double> potential;
  std::vector<double> force;
  std::vector<double> virial;
};

std::vector<int> parse_thread_list(const std::string& text)
{
  std::vector<int> threads;
  std::stringstream stream(text);
  std::string item;

  while (std::getline(stream, item, ',')) {
    const int value = std::atoi(item.c_str());
    if (value > 0) {
      threads.push_back(value);
    }
  }

  if (threads.empty()) {
    threads.push_back(1);
  }
  return threads;
}

void set_num_threads(const int threads)
{
#if defined(_OPENMP)
  omp_set_num_threads(threads);
#else
  (void)threads;
#endif
}

NepOutput compute_once(const std::string& model_file, const int threads)
{
  set_num_threads(threads);

  NEP nep;
  nep.init_from_file(model_file, false);

  const std::vector<int> type = {0, 1, 0, 1, 0, 1};
  const std::vector<double> box = {
    10.0, 0.0, 0.0,
    0.0, 10.0, 0.0,
    0.0, 0.0, 10.0};

  // position is stored as x[0..N), y[0..N), z[0..N), matching NEP::compute.
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

double max_abs_diff(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
  double result = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    result = std::max(result, std::abs(lhs[i] - rhs[i]));
  }
  return result;
}

double average_runtime_ms(const std::string& model_file, const int threads, const int repeats)
{
  const auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < repeats; ++i) {
    compute_once(model_file, threads);
  }
  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double, std::milli> elapsed = end - start;
  return elapsed.count() / repeats;
}

} // namespace

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " MODEL_FILE [THREADS=1,2,4] [REPEATS=20]\n";
    return 2;
  }

  const std::string model_file = argv[1];
  const std::vector<int> thread_list = parse_thread_list(argc > 2 ? argv[2] : "1,2,4");
  const int repeats = argc > 3 ? std::max(1, std::atoi(argv[3])) : 20;

  const NepOutput reference = compute_once(model_file, 1);

  std::cout << "threads,avg_ms,max_potential_diff,max_force_diff,max_virial_diff\n";
  for (const int threads : thread_list) {
    const NepOutput output = compute_once(model_file, threads);
    const double avg_ms = average_runtime_ms(model_file, threads, repeats);

    std::cout << threads << ","
              << avg_ms << ","
              << max_abs_diff(reference.potential, output.potential) << ","
              << max_abs_diff(reference.force, output.force) << ","
              << max_abs_diff(reference.virial, output.virial) << "\n";
  }

  return 0;
}
