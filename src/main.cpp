#include <iostream>
#include <vector>
#include <fftw3.h>
#include <sndfile.h>
#include <cmath>
#include <filesystem>

// Function to apply a Hamming window
void apply_hamming_window(std::vector<double>& window) {
    int N = window.size();
    for (int i = 0; i < N; ++i) {
        window[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (N - 1));
    }
}

// Function to compare FFT results
bool is_significantly_different(const fftw_complex* current, const fftw_complex* previous, int size, double threshold) {
    for (int i = 0; i < size; ++i) {
        double diff_real = std::abs(current[i][0] - previous[i][0]);
        double diff_imag = std::abs(current[i][1] - previous[i][1]);
        if (diff_real > threshold || diff_imag > threshold) {
            return true;
        }
    }
    return false;
}

// Function to find the dominant frequency
double find_dominant_frequency(const fftw_complex* out, int size, double sample_rate) {
    double max_magnitude = 0.0;
    int max_index = 0;
    for (int i = 0; i < size/2; ++i) {
        double magnitude = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            max_index = i;
        }
    }
    // Convert index to frequency
    double dominant_frequency = max_index * sample_rate / size;
    return dominant_frequency;
}

void process_audio_samples(const std::vector<double>& samples, int window_size, int hop_size, double sample_rate) {
    int N = samples.size();
    int num_windows = (N - window_size) / hop_size + 1;

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_complex *previous_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_plan plan = fftw_plan_dft_1d(window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    std::vector<double> window(window_size);
    bool first_window = true;
    double threshold = 0.1; // Set a threshold for significant difference

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

        // Compare with previous window and output if significantly different
        if (first_window || is_significantly_different(out, previous_out, window_size, threshold)) {
            double dominant_frequency = find_dominant_frequency(out, window_size, sample_rate);
            std::cout << "Window " << w << ": Dominant Frequency = " << dominant_frequency << " Hz\n";
            first_window = false;
        }

        // Copy current output to previous output
        for (int i = 0; i < window_size; ++i) {
            previous_out[i][0] = out[i][0];
            previous_out[i][1] = out[i][1];
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    fftw_free(previous_out);
}

int main() {
    // Open the audio file
    SF_INFO sfinfo;
    std::string file_path = "./tests/test2.wav";
    std::cout << "Attempting to open file: " << file_path << std::endl;
    SNDFILE *infile = sf_open(file_path.c_str(), SFM_READ, &sfinfo);
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
    double sample_rate = sfinfo.samplerate;
    // process_audio_samples(samples, window_size, hop_size, sample_rate);
    std::vector<double> mono_samples(sfinfo.frames);
    for (int i = 0; i < sfinfo.frames; ++i) {
        if (sfinfo.channels > 1) {
            double left = samples[i * sfinfo.channels];
            double right = samples[i * sfinfo.channels + 1];
            mono_samples[i] = (left + right) / 2.0;
        } else {
            mono_samples[i] = samples[i];
        }
    }
    // Use mono_samples instead of samples below
    process_audio_samples(mono_samples, window_size, hop_size, sample_rate);

    return 0;
}