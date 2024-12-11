#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

// Function declarations
void apply_hamming_window(std::vector<double> &data);
void fft_cooley_tukey(std::vector<std::complex<double>> &data);
void fft_cooley_tukey_parallel(std::vector<std::complex<double>> &data, int thread_count);

#endif // FFT_H