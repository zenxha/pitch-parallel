#include "fft.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <thread>

void apply_hamming_window(std::vector<double> &data) {
    size_t N = data.size();
    for (size_t i = 0; i < N; ++i) {
        data[i] *= 0.54 - 0.46 * std::cos(2 * M_PI * i / (N - 1));
    }
}

void fft_cooley_tukey(std::vector<std::complex<double>> &data) {
    size_t N = data.size();
    if (N <= 1) return;

    std::vector<std::complex<double>> even(N / 2), odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft_cooley_tukey(even);
    fft_cooley_tukey(odd);

    for (size_t k = 0; k < N / 2; ++k) {
        auto t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        data[k] = even[k] + t;
        data[k + N / 2] = even[k] - t;
    }
}

// Transforms time-domain signal to frequency-domain signal
void fft_cooley_tukey_parallel(std::vector<std::complex<double>> &data, int thread_count) {
  // Compute N and number of stages (log_2 (N))
  // Stages are not part of the DFT formula, but part of Cooley-Tukey optimization
  int N = data.size();
  int logN = std::ceil(std::log2(N));

  // Bit reversal
  for (int i = 1, j = 0; i < N; ++i) {
    int bit = N >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(data[i], data[j]);
    }
  }

  // Cooley-Tukey specific stage iteration
  for (int stage = 1; stage <= logN; ++stage) {
    // Compute group sizes; start with small 2-point DFT's, then combine to 4 point, etc.
    int group_size = 1 << stage;
    int group_size_half = group_size >> 1;

    // This value is equivalent to e^{-i * pi * (2^(stage-1))^(-1)} = e^{-2 * i * pi * (1/group_size)}
    // Or alternatively e^{-i * pi * (1/group_size_half)}
    // This is W(k=1,N=group_size)
    // Used later to compute W(k=j, group_size) = (base_factor)^j
    std::complex<double> base_factor = std::exp(std::complex<double>(0, -2.0 * M_PI / group_size));

    // Parallelism within each stage
    // Cannot parallelize multiple stages at once because data dependencies
    // Similar to hotplate assignment; Could not compute next iteration until dependencies were resolved
    auto group_process_lambda = [&](int start, int end) {
      // Compute the actual DFT from smaller DFTs
      // Smallest DFT (1 point) is just the input
      for (int i = start; i < end; i += group_size) {
        std::complex<double> rotating_factor = 1.0;  // W(k=0,N=group_size) = 1 (trivial)
        for (int j = 0; j < group_size_half; ++j) {
          std::complex<double> upper_term = rotating_factor * data[i + j + group_size_half];
          std::complex<double> lower_term = data[i + j];

          // Each lambda accesses data at indices:
          //  - start                   -> start + group_size_half - 1
          //  - start + group_size_half -> end - 1

          // intuitively, each iteration of the outer loop
          // for (int i = start; i < end; i += group_size)
          // touches from index i to index i + group_size
          // because this loop for (int j = 0; j < group_size_half; ++j)
          // accesses from i + 0 through i + 2*(group_size_half) - 1
          data[i + j] = lower_term + upper_term;
          data[i + j + group_size_half] = lower_term - upper_term;
          rotating_factor *= base_factor;
        }
      }
    };

    int group_count = N / group_size;
    int groups_per_thread = group_count / thread_count;
    int excess_groups = group_count % thread_count;

    std::vector<std::thread> threads;

    int start_group = 0;
    for (int thread_id = 0; thread_id < thread_count; ++thread_id) {
      int last_group = start_group + groups_per_thread;
      if (thread_id < excess_groups) {
        last_group += 1;
      }

      int group_start_index = start_group * group_size;
      int group_end_index = last_group * group_size;
      threads.emplace_back(group_process_lambda, group_start_index, group_end_index);

      start_group = last_group;
      // Causes next lambda to start at the index +1 the final one used in the previous lambda
    }

    for (auto &thread : threads) {
      thread.join();
    }
  }
}