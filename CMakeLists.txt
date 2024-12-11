cmake_minimum_required(VERSION 3.10)
project(pitch-parallel)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Define the executable target
add_executable(main src/main.cpp src/fft.cpp)

# Link libraries
if(OpenMP_CXX_FOUND)
    target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
endif()

# Include FFTW, sndfile for macOS
target_include_directories(main PRIVATE /opt/homebrew/Cellar/fftw/3.3.10_1/include)
target_link_directories(main PRIVATE /opt/homebrew/Cellar/fftw/3.3.10_1/lib)

target_include_directories(main PRIVATE /opt/homebrew/Cellar/libsndfile/1.2.2/include)
target_link_directories(main PRIVATE /opt/homebrew/Cellar/libsndfile/1.2.2/lib)

target_include_directories(main PRIVATE /opt/homebrew/Cellar/portaudio/19.7.0/include)
target_link_directories(main PRIVATE /opt/homebrew/Cellar/portaudio/19.7.0/lib)

target_link_libraries(main PRIVATE fftw3 sndfile portaudio)
