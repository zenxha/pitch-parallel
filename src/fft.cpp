#include "fft.h"
#include <cmath>
#include <algorithm>
#include <thread>
#define M_PI 3.14159265358979323846264338327950288

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

double compute_threshold_dynamic(const std::vector<std::complex<double>> &fft_result, double factor) {
    double sum = 0.0;
    for (const auto &val : fft_result) {
        sum += std::abs(val);
    }
    return factor * (sum / fft_result.size());
}

std::vector<std::pair<double, double>> find_prominent_frequencies(
    const std::vector<std::complex<double>> &fft_result, double sample_rate, double threshold) {
    std::vector<std::pair<double, double>> frequencies;
    size_t N = fft_result.size();

    for (size_t i = 0; i < N / 2; ++i) {
        double magnitude = std::abs(fft_result[i]);
        if (magnitude > threshold) {
            double frequency = i * sample_rate / N;
            frequencies.emplace_back(frequency, magnitude);
        }
    }

    std::sort(frequencies.begin(), frequencies.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.second > b.second; });

    return frequencies;
}

double find_dominant_frequency(const std::vector<std::complex<double>> &spectrum, double sample_rate) {
  int N = (int)spectrum.size();
  // Only consider first N/2 (real spectrum)
  int half = N / 2;
  double max_magnitude = 0.0;
  int max_index = 0;
  for (int i = 0; i < half; ++i) {
    double mag = std::abs(spectrum[i]);
    if (mag > max_magnitude) {
      max_magnitude = mag;
      max_index = i;
    }
  }
  return (max_index * sample_rate) / N;
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


std::string frequencyToPitchName(double frequency, double referenceFrequency) {
    const std::vector<std::string> noteNames = {"C ", "C#", "D ", "D#", "E ", "F ", "F#", "G ", "G#", "A ", "A#", "B "};
    
    // Calculate the number of semitones relative to A4
    double pitch = 12 * std::log2(frequency / referenceFrequency);
    
    // Find the nearest MIDI note number
    int midiNote = static_cast<int>(std::round(pitch)) + 69;
    
    // Determine the octave (MIDI note 0 is C-1)
    int octave = (midiNote / 12) - 1;
    
    // Determine the note name
    std::string noteName = noteNames[midiNote % 12];
    
    // Combine note name and octave
    return noteName + std::to_string(octave);
}