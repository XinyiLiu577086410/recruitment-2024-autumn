#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "SmithWaterman.hpp"
#include "Timer.hpp"

int main(int argc, char **argv) {
  ScopeTimer timer("Elapsed time");

  if (argc != 4) {
    std::cerr << "Number of parameters is not correct" << std::endl;
    std::cerr << "Usage: " << argv[0]
              << " query_seq_path target_seq_path ref_path" << std::endl;
    throw std::runtime_error("Incorrect number of parameters");
  }

  std::filesystem::path query_seq_path = argv[1];
  std::filesystem::path target_seq_path = argv[2];
  std::filesystem::path ref_path = argv[3];

  SmithWaterman solver = SmithWaterman(query_seq_path, target_seq_path);

  solver.solve();

  auto res = solver.validate(ref_path);

  return res;
}
