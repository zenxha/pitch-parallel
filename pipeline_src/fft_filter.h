#ifndef __FFT_FILTER_H__
#define __FFT_FILTER_H__

#include <vector>
#include <complex>

void apply_band_pass_filter(std::vector<std::complex<double>> &spectrum, double sample_rate, double low_cutoff, double high_cutoff);

#endif // __FFT_FILTER_H__