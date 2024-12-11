#ifndef __FFT_ANALYSIS_H__
#define __FFT_ANALYSIS_H__

#include <algorithm>
#include <complex>
#include <string>
#include <vector>

double compute_threshold_dynamic(const std::vector<std::complex<double>> &fft_result, double factor);
std::vector<std::pair<double, double>> find_prominent_frequencies(const std::vector<std::complex<double>> &fft_result, double sample_rate, double threshold);
double find_dominant_frequency(const std::vector<std::complex<double>> &spectrum, double sample_rate);
// Function to convert a frequency to a pitch name
std::string frequencyToPitchName(double frequency, double referenceFrequency = 440.0);

#endif