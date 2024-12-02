#include <iostream>
#include <omp.h>
#include <fftw3.h>


int main() {
    // Define the size of the transform
    const int N = 8;

    // Allocate input and output arrays using fftw_malloc
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Create a plan for the forward FFT
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Initialize input data (hard-coded example)
    for (int i = 0; i < N; ++i) {
        in[i][0] = i + 1; // Real part
        in[i][1] = 0.0;   // Imaginary part
    }

    // Execute the FFT
    fftw_execute(plan);

    // Output the results
    std::cout << "Input:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "in[" << i << "] = " << in[i][0] << " + " << in[i][1] << "i" << std::endl;
    }

    std::cout << "\nOutput (FFT):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "out[" << i << "] = " << out[i][0] << " + " << out[i][1] << "i" << std::endl;
    }

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return 0;
}
