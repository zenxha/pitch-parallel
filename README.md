# pitch-parallel

This project uses a the Cooley-Tukey Fourier Transform (FFT) with multithreading to identify musical notes in real-time and prerecorded audio samples by analyzing frequency components. 

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
