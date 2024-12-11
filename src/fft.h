#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

// Function declarations
void apply_hamming_window(std::vector<double> &data);
void fft_cooley_tukey(std::vector<std::complex<double>> &data);
void fft_cooley_tukey_parallel(std::vector<std::complex<double>> &data, int thread_count);
double compute_threshold_dynamic(const std::vector<std::complex<double>> &fft_result, double factor);
std::vector<std::pair<double, double>> find_prominent_frequencies(
const std::vector<std::complex<double>> &fft_result, double sample_rate, double threshold);
double find_dominant_frequency(const std::vector<std::complex<double>> &spectrum, double sample_rate);
// Function to convert a frequency to a pitch name
std::string frequencyToPitchName(double frequency, double referenceFrequency = 440.0);


#endif // FFT_H