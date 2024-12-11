#include <sndfile.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <thread>
#include <vector>

#include "fftw3.h"

void apply_hamming_window(std::vector<double> &window) {
  int N = window.size();
  for (int i = 0; i < N; ++i) {
    window[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (N - 1));
  }
}

double compute_threshold(const std::vector<std::complex<double>> &spectrum, double factor = 2.0) {
  int N = spectrum.size() / 2;  // Use only the positive frequencies
  std::vector<double> magnitudes(N);
  for (int i = 0; i < N; ++i) {
    magnitudes[i] = std::abs(spectrum[i]);
  }
  std::nth_element(magnitudes.begin(), magnitudes.begin() + N / 2, magnitudes.end());
  double median_magnitude = magnitudes[N / 2];
  return factor * median_magnitude;
}

double compute_threshold_dynamic(const std::vector<std::complex<double>> &spectrum, double percentage = 0.1) {
  int N = spectrum.size() / 2;  // Use only the positive frequencies
  double max_magnitude = 0.0;
  for (int i = 0; i < N; ++i) {
    max_magnitude = std::max(max_magnitude, std::abs(spectrum[i]));
  }
  return percentage * max_magnitude;
}

void fft_cooley_tukey(std::vector<std::complex<double>> &data) {
  int N = (int)data.size();
  int logN = 0;
  while ((1 << logN) < N) {
    logN++;
  }

  // Bit-reversal
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

  // Iterative FFT
  for (int s = 1; s <= logN; ++s) {
    int m = 1 << s;
    int m2 = m >> 1;
    std::complex<double> w_m = std::exp(std::complex<double>(0, -M_PI / m2));
    for (int k = 0; k < N; k += m) {
      std::complex<double> w = 1.0;
      for (int j = 0; j < m2; ++j) {
        std::complex<double> t = w * data[k + j + m2];
        std::complex<double> u = data[k + j];
        data[k + j] = u + t;
        data[k + j + m2] = u - t;
        w *= w_m;
      }
    }
  }
}

// Function to find the dominant frequency
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
// Suppose these are already defined somewhere:
// void fft_cooley_tukey(std::vector<std::complex<double>> &data);
// double find_dominant_frequency(const std::vector<std::complex<double>> &spectrum, double sample_rate);
// void apply_hamming_window(std::vector<double> &window);

// Finds prominent frequencies in the spectrum
std::vector<std::pair<double, double>> find_prominent_frequencies(const std::vector<std::complex<double>> &spectrum, double sample_rate, double threshold) {
  int N = spectrum.size();
  int half = N / 2;  // Only consider first half (positive frequencies)

  std::vector<std::pair<double, double>> prominent_frequencies;

  for (int i = 1; i < half - 1; ++i) {
    // Compute magnitude
    double magnitude = std::abs(spectrum[i]);

    // std::cout << magnitude << std::endl;

    // Check if this is a local maximum and above the threshold
    if (magnitude > threshold && magnitude > std::abs(spectrum[i - 1]) && magnitude > std::abs(spectrum[i + 1])) {
      double frequency = (i * sample_rate) / N;  // Convert bin index to frequency
      prominent_frequencies.emplace_back(frequency, magnitude);
    }
  }

  // Sort by magnitude (descending order)
  std::sort(prominent_frequencies.begin(), prominent_frequencies.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b) { return a.second > b.second; });

  return prominent_frequencies;
}

void process_window_range(int start_window, int end_window, const std::vector<double> &mono, int window_size, int hop_size, double sample_rate, std::vector<std::pair<int, std::vector<std::pair<double, double>>>> &thread_results, double threshold) {
  for (int w = start_window; w < end_window; ++w) {
    // Extract current window
    std::vector<double> window(window_size);
    for (int i = 0; i < window_size; ++i) {
      window[i] = mono[w * hop_size + i];
    }

    // Apply Hamming window
    apply_hamming_window(window);

    // Prepare complex vector for FFT
    std::vector<std::complex<double>> fft_input(window_size);
    for (int i = 0; i < window_size; ++i) {
      fft_input[i] = std::complex<double>(window[i], 0.0);
    }

    // Compute FFT
    fft_cooley_tukey(fft_input);

    // Find prominent frequencies
    auto computed_threshold = compute_threshold_dynamic(fft_input, 0.35);
    // std::cout << w << " : Selected threshold = " << computed_threshold << std::endl;
    auto prominent_frequencies = find_prominent_frequencies(fft_input, sample_rate, computed_threshold);

    // Store result in thread-local vector
    thread_results.emplace_back(w, prominent_frequencies);
  }
}
void parallel_process_audio(const std::vector<double> &mono, int window_size, int hop_size, double sample_rate, int num_windows, double threshold) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) num_threads = 2;  // fallback if hardware_concurrency is unreliable

  int windows_per_thread = num_windows / num_threads;
  int remainder = num_windows % num_threads;

  std::vector<std::thread> threads;
  std::vector<std::vector<std::pair<int, std::vector<std::pair<double, double>>>>> thread_results(num_threads);
  // std::vector<std::vector<std::pair<int, std::pair<double, double>>>> thread_results(num_threads);

  int current_start = 0;

  // Launch threads
  for (unsigned int t = 0; t < num_threads; ++t) {
    int current_end = current_start + windows_per_thread + (t < remainder ? 1 : 0);
    threads.emplace_back(process_window_range, current_start, current_end, std::ref(mono), window_size, hop_size, sample_rate, std::ref(thread_results[t]), threshold);
    current_start = current_end;
  }

  // Join threads
  for (auto &thr : threads) {
    thr.join();
  }

  // Merge results from all threads
  std::vector<std::pair<int, std::vector<std::pair<double, double>>>> results;
  for (const auto &thread_result : thread_results) {
    results.insert(results.end(), thread_result.begin(), thread_result.end());
  }

  // Sort results by window index
  std::sort(results.begin(), results.end());

  // Output sorted results
  for (const auto &result : results) {
    // double timestamp = static_cast<double>(result.first * hop_size) / sample_rate;
    double window_center_offset = static_cast<double>(window_size / 2) / sample_rate;
    double timestamp = static_cast<double>(result.first * hop_size) / sample_rate + window_center_offset;

    std::cout << "Window " << result.first << " (at " << timestamp << " s): ";
    // std::cout << "Window " << result.first << ": ";
    // std::cout << result.second.size () << " prominent frequencies: ";
    for (const auto &freq : result.second) {
      std::cout << freq.first << " Hz (" << freq.second << "), ";
    }
    std::cout << "\n";
  }
}

int main() {
  // Open the audio file using libsndfile
  SF_INFO sfinfo;
  std::string file_path = "./tests/test1.wav";
  std::cout << "Attempting to open file: " << file_path << std::endl;
  SNDFILE *infile = sf_open(file_path.c_str(), SFM_READ, &sfinfo);
  if (!infile) {
    std::cerr << "Error opening audio file!" << std::endl;
    return 1;
  }

  // Read stereo samples
  std::vector<double> samples(sfinfo.frames * sfinfo.channels);
  sf_readf_double(infile, samples.data(), sfinfo.frames);
  sf_close(infile);

  // Convert to mono by averaging channels
  std::vector<double> mono(sfinfo.frames);
  for (int i = 0; i < sfinfo.frames; ++i) {
    if (sfinfo.channels == 2) {
      double left = samples[i * sfinfo.channels];
      double right = samples[i * sfinfo.channels + 1];
      mono[i] = (left + right) / 2.0;
    } else {
      mono[i] = samples[i];
    }
  }

  // Parameters for short-time analysis
  int window_size = 4096;  // Must be power of two for radix-2 FFT
  int hop_size = 512;
  double sample_rate = sfinfo.samplerate;
  int N = (int)mono.size();
  int num_windows = (N - window_size) / hop_size + 1;
  std::cout << "N = " << N << std::endl;
  std::cout << "num_windows = " << num_windows << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  parallel_process_audio(mono, window_size, hop_size, sample_rate, num_windows, 10.0);

  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout << "Processing took " << duration << " milliseconds." << std::endl;

  return 0;
}