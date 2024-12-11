#include "fft_filter.h"

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