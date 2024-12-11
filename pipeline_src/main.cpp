#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include "fft.h"
#include "fft_analysis.h"
#include "fft_filter.h"
#include "fft_format.h"
#include "portaudio.h"
#include "sndfile.h"

// Constants / Parameters
constexpr int SAMPLE_RATE = 44100;
constexpr int FRAMES_PER_BUFFER = 512;
constexpr int NUM_CHANNELS = 1;
int FFT_WINDOW_SIZE = 1024;  // power of 2 is required for radix-2 FFT
constexpr int HOP_SIZE = 512;
constexpr int LOW_CUTOFF = 20;
constexpr int HIGH_CUTOFF = 5000;

// Control Atomics
std::atomic<bool> exit_flag(false);  // for exiting the main thread

// Types / Classes
typedef std::complex<double> cDouble;
typedef std::vector<cDouble> FFTData;

struct AudioData {
  std::vector<double> buffer;
  std::mutex mutex;
};

template <typename T>
class ThreadSafeQueue {
 private:
  std::queue<T> q;
  std::mutex mtx;
  std::condition_variable cv;

 public:
  void push(T item) {
    std::lock_guard<std::mutex> lock(mtx);
    q.push(std::move(item));
    cv.notify_one();
  }

  T pop() {
    std::unique_lock<std::mutex> lock(mtx);
    cv.wait(lock, [this] { return !q.empty(); });
    T item = std::move(q.front());
    q.pop();
    return item;
  }
};

// Portaudio Helpers
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

// Thread Functions (Task Level Parallelism)
// Handles moving data from audio buffer (buffered mic input)
// into the TODO queue as a window.
void framing_thread_func(AudioData &audio_data, ThreadSafeQueue<std::vector<double>> &time_domain_queue, int FFT_WINDOW_SIZE, int HOP_SIZE) {
  while (!exit_flag) {
    std::vector<double> window;

    {
      //
      std::lock_guard<std::mutex> lock(audio_data.mutex);
      if (audio_data.buffer.size() >= FFT_WINDOW_SIZE) {
        // Move audio data from buffer to window
        window.assign(audio_data.buffer.begin(), audio_data.buffer.begin() + FFT_WINDOW_SIZE);
        audio_data.buffer.erase(audio_data.buffer.begin(), audio_data.buffer.begin() + HOP_SIZE);
      }
    }

    if (!window.empty()) {
      apply_hamming_window(window);
      time_domain_queue.push(std::move(window));
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
  }
}

// Thread that fetches work to perform FFT on
void fft_thread_func(ThreadSafeQueue<std::vector<double>> &time_domain_queue, ThreadSafeQueue<FFTData> &freq_domain_queue, bool use_parallel, int num_threads, double sample_rate) {
  while (!exit_flag) {
    std::vector<double> window = time_domain_queue.pop();

    FFTData fft_input(window.size());
    for (size_t i = 0; i < window.size(); ++i) {
      fft_input[i] = cDouble(window[i], 0.0);
    }

    if (use_parallel) {
      fft_cooley_tukey_parallel(fft_input, num_threads);
    } else {
      fft_cooley_tukey(fft_input);
    }

    apply_band_pass_filter(fft_input, sample_rate, LOW_CUTOFF, HIGH_CUTOFF);
    freq_domain_queue.push(std::move(fft_input));
  }
}

void analysis_thread_func(ThreadSafeQueue<FFTData> &freq_domain_queue, double sample_rate) {
  std::string last_output;
  while (!exit_flag) {
    auto fft_output = freq_domain_queue.pop();

    auto threshold = compute_threshold_dynamic(fft_output, 0.35);
    auto dominant_frequency = find_dominant_frequency(fft_output, sample_rate);
    auto pitch = frequencyToPitchName(dominant_frequency);

    // Print results
    std::ostringstream oss;
    oss << pitch << "\t" << formatDouble7(dominant_frequency) << " Hz";
    std::string output = oss.str();
    std::cout << "\r" << std::string(last_output.size(), ' ') << "\r";
    std::cout << output << std::flush;
    last_output = output;
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

  ThreadSafeQueue<std::vector<double>> time_domain_queue;
  ThreadSafeQueue<std::vector<std::complex<double>>> freq_domain_queue;

  std::cout << "Should we use parallel? (y?)" << std::endl;

  char input;
  std::cin >> input;

  bool use_parallel = false;
  if (input == 'y') {
    use_parallel = true;
  }

  int num_threads = 1;
  if (use_parallel) {
    std::cout << "\nHow many threads do you want to use? (Max 16, min 1)" << std::endl;
    std::cin >> num_threads;

    if (num_threads > 16) {
      return -1;
    } else if (num_threads < 1) {
      return -1;
    }
  }

  std::cout << "\nWindow Size = 1024 * K" << std::endl << "K = ";
  int K;
  std::cin >> K;
  if (K < 1) {
    return -1;
  } else if (K > 16) {
    return -1;
  }
  FFT_WINDOW_SIZE = 1024 * K;

  std::thread frame_thread(framing_thread_func, std::ref(audio_data), std::ref(time_domain_queue), FFT_WINDOW_SIZE, HOP_SIZE);
  std::thread fft_thread(fft_thread_func, std::ref(time_domain_queue), std::ref(freq_domain_queue), use_parallel, num_threads, SAMPLE_RATE);
  std::thread analysis_thread(analysis_thread_func, std::ref(freq_domain_queue), SAMPLE_RATE);

  std::cout << "\nPress 'q + Enter' to quit." << std::endl;
  char exit_input;
  while (std::cin >> exit_input) {
    if (exit_input == 'q') {
      exit_flag = true;  // Signal the threads to stop
      break;
    }
  }

  frame_thread.join();
  fft_thread.join();
  analysis_thread.join();

  Pa_StopStream(stream);
  Pa_CloseStream(stream);
  Pa_Terminate();

  return 0;
}