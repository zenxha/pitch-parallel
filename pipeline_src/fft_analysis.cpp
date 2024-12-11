#include "fft_analysis.h"

double compute_threshold_dynamic(const std::vector<std::complex<double>> &fft_result, double factor) {
  double sum = 0.0;
  for (const auto &val : fft_result) {
    sum += std::abs(val);
  }
  return factor * (sum / fft_result.size());
}

std::vector<std::pair<double, double>> find_prominent_frequencies(const std::vector<std::complex<double>> &fft_result, double sample_rate, double threshold) {
  std::vector<std::pair<double, double>> frequencies;
  size_t N = fft_result.size();

  for (size_t i = 0; i < N / 2; ++i) {
    double magnitude = std::abs(fft_result[i]);
    if (magnitude > threshold) {
      double frequency = i * sample_rate / N;
      frequencies.emplace_back(frequency, magnitude);
    }
  }

  std::sort(frequencies.begin(), frequencies.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b) { return a.second > b.second; });

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