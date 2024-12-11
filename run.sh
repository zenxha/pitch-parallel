#!/bin/bash

# Path to the audio file
FILE_PATH="./tests/piano.wav"

# Number of tests to run for each thread count
NUM_TESTS=5
OUTPUT_FILE="results.txt"

# Clear the output file if it exists
> $OUTPUT_FILE

# Loop over thread counts from 1 to 8
for THREAD_COUNT in {1..8}; do
    echo "Running tests with $THREAD_COUNT threads..." | tee -a $OUTPUT_FILE
    # Run the specified number of tests for each thread count
    for TEST_NUM in $(seq 1 $NUM_TESTS); do
        echo "Test $TEST_NUM with $THREAD_COUNT threads:" | tee -a $OUTPUT_FILE
        # Capture the output time and append it to the output file
        TIME_OUTPUT=$(./build/main2 $THREAD_COUNT $FILE_PATH)
        echo $TIME_OUTPUT | tee -a $OUTPUT_FILE
    done
    echo "" | tee -a $OUTPUT_FILE
done

echo "All tests completed." | tee -a $OUTPUT_FILE