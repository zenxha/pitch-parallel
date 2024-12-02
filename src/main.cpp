#include <iostream>
#include <vector>
#include <fftw3.h>
#include <sndfile.h>
#include <cmath>

// Function to apply a Hamming window
void apply_hamming_window(std::vector<double>& window) {
    int N = window.size();
    for (int i = 0; i < N; ++i) {
        window[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (N - 1));
    }
}

void process_audio_samples(const std::vector<double>& samples, int window_size, int hop_size) {
    int N = samples.size();
    int num_windows = (N - window_size) / hop_size + 1;

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_plan plan = fftw_plan_dft_1d(window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    std::vector<double> window(window_size);

    for (int w = 0; w < num_windows; ++w) {
        // Copy samples into window and apply Hamming window
        for (int i = 0; i < window_size; ++i) {
            window[i] = samples[w * hop_size + i];
        }
        apply_hamming_window(window);

        // Prepare input for FFT
        for (int i = 0; i < window_size; ++i) {
            in[i][0] = window[i];
            in[i][1] = 0.0;
        }

        // Execute FFT
        fftw_execute(plan);

        // Output the results for this window
        std::cout << "Window " << w << ":\n";
        for (int i = 0; i < window_size; ++i) {
            std::cout << "out[" << i << "] = " << out[i][0] << " + " << out[i][1] << "i\n";
        }
        std::cout << "\n";
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

int main() {
    // Open the audio file
    SF_INFO sfinfo;
    SNDFILE *infile = sf_open("./tests/test1.wav", SFM_READ, &sfinfo);
    if (!infile) {
        std::cerr << "Error opening audio file!" << std::endl;
        return 1;
    }

    // Read the audio samples
    std::vector<double> samples(sfinfo.frames * sfinfo.channels);
    sf_readf_double(infile, samples.data(), sfinfo.frames);
    sf_close(infile);

    // Process the audio samples
    int window_size = 1024;
    int hop_size = 512;
    process_audio_samples(samples, window_size, hop_size);

    return 0;
}