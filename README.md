# pitch-parallel

This project uses a Discrete Fourier Transform (DFT) with OpenMP to identify musical notes in audio samples by analyzing frequency components in parallel. Ideal for musicians and music learners, this tool isolates dominant frequencies in real-time and maps them to musical notes for easy song analysis.

## Environment Setup Instructions for Mac

1. **Install Dependencies:**
    ```sh
    brew install fftw libomp
    ```

2. **Set Environment Variables:**
    ```sh
    export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
    export LDFLAGS="-L/opt/homebrew/Cellar/fftw/3.3.10_1/lib"
    export CPPFLAGS="-I/opt/homebrew/Cellar/fftw/3.3.10_1/include"
    ```

3. **Compile the Project:**
    ```sh
    g++ -I/usr/local/include -L/usr/local/lib -lfftw3 -o main <main.cpp path>
    ```

4. **Build with CMake:**
    ```sh
    cmake -S . -B build \
         -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
         -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
    cmake --build build
    ```
