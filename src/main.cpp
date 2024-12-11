#include <algorithm>
#include <chrono>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>
#include <thread>
#include <vector>
#include <iomanip>
#include <atomic>
#include <mutex>

#include "fft.h"
#include "portaudio.h"
#include "sndfile.h"

constexpr int SAMPLE_RATE = 44100;
constexpr int FRAMES_PER_BUFFER = 512;
constexpr int NUM_CHANNELS = 1;
constexpr int HOP_SIZE = 512;
constexpr int LOW_CUTOFF = 20;
constexpr int HIGH_CUTOFF = 5000;
int FFT_WINDOW_SIZE = 1024;  // power of 2 is required for radix-2 FFT

std::atomic<bool> exit_flag(false); // for exiting the main thread

struct AudioData {
  std::vector<double> buffer;
  std::mutex mutex;
};

int audio_callback(const void *input_buffer, void *output_buffer, unsigned long frames_per_buffer, const PaStreamCallbackTimeInfo *time_info, PaStreamCallbackFlags status_flags, void *user_data) {
  auto *audio_data = static_cast<AudioData *>(user_data);

  // copy audio input to buffer
  const float *input = static_cast<const float *>(input_buffer);
  std::lock_guard<std::mutex> lock(audio_data->mutex);
  for (unsigned long i = 0; i < frames_per_buffer; ++i) {
    audio_data->buffer.push_back(input[i]);
  }

  return paContinue;
}

void apply_band_pass_filter(std::vector<std::complex<double>> &spectrum, double sample_rate, double low_cutoff, double high_cutoff) {
  int N = spectrum.size();
  int low_bin = (low_cutoff / sample_rate) * N;
  int high_bin = (high_cutoff / sample_rate) * N;

  for (int i = 0; i < N; ++i) {
    if (i < low_bin || i > high_bin) {
      spectrum[i] = 0;
    }
  }
}

std::string formatDouble7(double value) {
    std::ostringstream oss;

    if (value >= 100000 || value <= -10000) {
        return " ERROR ";
    }

    if (std::abs(value) >= 1000) {
        oss << std::fixed << std::setprecision(2);
    } else if (std::abs(value) >= 100) {
        oss << std::fixed << std::setprecision(3);
    } else if (std::abs(value) >= 10) {
        oss << std::fixed << std::setprecision(4);
    } else {
        oss << std::fixed << std::setprecision(5);
    }

    oss << value;

    std::string result = oss.str();

    if (result.length() > 7) {
        result = result.substr(0, 7);
    } else if (result.length() < 7) {
        result.insert(0, 7 - result.length(), ' ');
    }

    return result;
}

void process_realtime_audio(AudioData &audio_data, double sample_rate, bool use_parallel, int num_threads) {
  std::string last_output;
  while (!exit_flag) {
    std::vector<double> window;

    {
      std::lock_guard<std::mutex> lock(audio_data.mutex);
      if (audio_data.buffer.size() >= FFT_WINDOW_SIZE) {
        window.assign(audio_data.buffer.begin(), audio_data.buffer.begin() + FFT_WINDOW_SIZE);
        audio_data.buffer.erase(audio_data.buffer.begin(), audio_data.buffer.begin() + HOP_SIZE);
      }
    }

    if (window.empty()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      continue;
    }

    auto start = std::chrono::high_resolution_clock::now();

    apply_hamming_window(window);

    std::vector<std::complex<double>> fft_input(FFT_WINDOW_SIZE);
    for (size_t i = 0; i < window.size(); ++i) {
      fft_input[i] = std::complex<double>(window[i], 0.0);
    }

   
    if (use_parallel) {
      fft_cooley_tukey_parallel(fft_input, num_threads);
    } else {
      fft_cooley_tukey(fft_input);
    }

    apply_band_pass_filter(fft_input, sample_rate, LOW_CUTOFF, HIGH_CUTOFF);

    // auto threshold = compute_threshold_dynamic(fft_input, 0.35);
    // auto prominent_frequencies = find_prominent_frequencies(fft_input, sample_rate, threshold);
    auto dominant_frequency = find_dominant_frequency(fft_input, SAMPLE_RATE);
    auto pitch = frequencyToPitchName(dominant_frequency);

    // std::cout << "Prominent frequencies:\n";
    // for (const auto &freq : prominent_frequencies) {
    //     std::cout << freq.first << " Hz (" << freq.second << "), ";
    // }
    auto end = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // One-line output
    std::ostringstream oss;
    oss << "[" << duration.count() << " ms]\t" << pitch << "\t" << formatDouble7(dominant_frequency) << " Hz";
    std::string output = oss.str();
    std::cout << "\r" << std::string(last_output.size(), ' ') << "\r";
    std::cout << output << std::flush;
    last_output = output;

    // Continuous output
    // std::cout << "[" << duration.count() << " ms] " << pitch << "\t" << formatDouble7(dominant_frequency) << " Hz" << std::endl;
  }
  std::cout << std::endl;
}

int main() {
  Pa_Initialize();

  AudioData audio_data;
  audio_data.buffer.reserve(FFT_WINDOW_SIZE * 2);  // buffer size to store overlapping windows

  PaStream *stream;
  Pa_OpenDefaultStream(&stream, NUM_CHANNELS, 0, paFloat32, SAMPLE_RATE, FRAMES_PER_BUFFER, audio_callback, &audio_data);
  Pa_StartStream(stream);

  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!Welcome to Pitch Parallel!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
  std::cout << "[!] Do you want to use parallel processing? (y/n): ";

  char input;
  std::cin >> input;

  bool use_parallel = input == 'y';
  int num_threads = 1;
  if (use_parallel) {
    std::cout << "\n[!] How many threads do you want to use? (max 16, min 1): ";
    std::cin >> num_threads;
    if (num_threads > 16 || num_threads < 1) {
      std::cout << "\nError: # of threads is out of range!" << std::endl;
      return -1;
    }
  }

  std::cout << "\n[?] The FFT window size is determined by the formula: 1024 * K";
  std::cout << "\n[!] Please enter a value for K (must be a power of 2): ";
  
  int K;
  std::cin >> K;
  if (K < 1 || K > 16) {
    std::cout << "\nError: K is out of range!" << std::endl;
    return -1;
  } else if ((K & (K - 1)) != 0){
    std::cout << "\nError: K must be a power of 2!" << std::endl;
    return -1;
  }
  FFT_WINDOW_SIZE = 1024 * K;

  std::thread processing_thread(process_realtime_audio, std::ref(audio_data), SAMPLE_RATE, use_parallel, num_threads);
  // process_realtime_audio(audio_data, SAMPLE_RATE, use_parallel, num_threads);

  std::cout << "\n[!] Press 'q + Enter' to quit.\n" << std::endl;
  std::cout << "Latency\tPitch\tFrequency" << std::endl;
  char exit_input;
  while (std::cin >> exit_input) {
    if (exit_input == 'q') {
      exit_flag = true;  // Signal the threads to stop
      break;
    }
  }

  processing_thread.join();
  std::cout << "\n[!] Quit." << std::endl;

  Pa_StopStream(stream);
  Pa_CloseStream(stream);
  Pa_Terminate();

  return 0;
}