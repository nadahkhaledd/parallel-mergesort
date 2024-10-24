#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <sys/time.h>
#include <iostream>
#include <algorithm>

#include <omp.h>
#include <cstdlib>
#include <cstdio>

#include <cmath>
#include <ctime>
#include <cstring>
#include <cctype>

using namespace std;

// Machine-specific settings
#define CACHE_LINE 64        
#define MIN_PARALLEL_SIZE 1024   
#define MAX_DEPTH 5 // 64 tasks  

/*
    @author nada-khaled
*/

/**
 * Verifies if the array is correctly sorted
 * @param ref Original array (sorted for comparison)
 * @param data Array to verify
 * @param size Length of the arrays
 * @return true if correctly sorted, false otherwise
 */
bool isSorted(int ref[], int data[], const size_t size) {
    sort(ref, ref + size);
    for (size_t idx = 0; idx < size; ++idx) {
        if (ref[idx] != data[idx]) {
            return false;
        }
    }
    return true;
}

/**
 * sequential merge step (straight-forward implementation)
 * @param out Output array where merged result will be stored
 * @param in Input array with two sorted sections
 * @param begin1 Start index of first sorted section
 * @param end1 End index of first sorted section
 * @param begin2 Start index of second sorted section
 * @param end2 End index of second sorted section
 * @param outBegin Starting position in output merged array
 */
void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
    long left = begin1;    
    long right = begin2;   
    long idx = outBegin;   

    while (left < end1 && right < end2) {
        if (in[left] <= in[right]) {
            out[idx] = in[left];
            left++;
        } else {
            out[idx] = in[right];
            right++;
        }
        idx++;
    }

    while (left < end1) {
        out[idx] = in[left];
        left++, idx++;
    }

    while (right < end2) {
        out[idx] = in[right];
        right++, idx++;
    }
}

/**
 * Sequential implementation of MergeSort
 * @param array Array to be sorted
 * @param tmp Temporary array 
 * @param inplace to sort array inplace
 * @param begin Start index of section to sort
 * @param end End index of section to sort
 */
void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
    if (begin < (end - 1)) {
        const long half = (begin + end) / 2;  
        MsSequential(array, tmp, !inplace, begin, half);
        MsSequential(array, tmp, !inplace, half, end);
        
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }
    } else if (!inplace) {
        tmp[begin] = array[begin]; 
    }
}

/**
 * Parallel implementation of MergeSort using OpenMP 
 * @param array Array to be sorted
 * @param tmp Temporary array
 * @param inplace to sort array inplace
 * @param begin Start index of section to sort
 * @param end End index of section to sort
 * @param depth Current recursion depth for controlling parallel tasks
 */
void MsMergeParallel(int *array, int *tmp, bool inplace, long begin, long end, int depth) {
    const long size = end - begin;

	/*#pragma omp critical
    {
        cout << "Creating task at depth " << depth << " for size " << size << endl;
    }*/
    
    // Switch to sequential sort for small arrays
    if (size < MIN_PARALLEL_SIZE/(1 << depth) || depth >= MAX_DEPTH) {
        if (inplace) {
            sort(array + begin, array + end);
        } else {
            sort(array + begin, array + end);
            memcpy(tmp + begin, array + begin, size * sizeof(int));
        }
        return;
    }

    if (begin < (end - 1)) {
        const long half = (begin + end) / 2;
        
        // Create parallel tasks for each half (Divide)
        #pragma omp task shared(array, tmp) if(size > MIN_PARALLEL_SIZE/(1 << depth))
        MsMergeParallel(array, tmp, !inplace, begin, half, depth + 1);
        
        #pragma omp task shared(array, tmp) if(size > MIN_PARALLEL_SIZE/(1 << depth))
        MsMergeParallel(array, tmp, !inplace, half, end, depth + 1);
        
        // wait to merge tasks
        #pragma omp taskwait

        // Merge final parts (Conquer)
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }
    }
}

/**
 * parallel merge sort execution
 * @param array Array to be sorted
 * @param tmp Temporary array
 * @param size Length of the array
 */
void MsParallel(int *array, int *tmp, const size_t size, int num_threads) {
    //Adjust number of threads based on array size
    //int num_threads = omp_get_max_threads();

    if (size < MIN_PARALLEL_SIZE * 4) {
        num_threads = std::min(num_threads, 2);
    }

    cout << "Start parallel implementation with  " << num_threads << " threads.\n";
    // Create parallel region with adjusted thread number
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp single // to ensure calling once by first thread
        {
            MsMergeParallel(array, tmp, true, 0, size, 0); // depth=0 as first step to increment one by one
        }
    }
}

/**
 * Calculates speedup between sequential and  parallel implementation
 * @param sequentialTime Sequential sort time
 * @param parallelTime Parallel sort time
 */
void calculateSpeedup(double sequentialTime, double parallelTime) {
    if (parallelTime > 0) {
        double speedup = sequentialTime / parallelTime;
        printf("Speedup: %f\n", speedup);
    } else {
        printf("Speedup calculation error: Parallel time = 0.\n");
    }
}

/**
 * check for sequential mode choices from terminal
 * @param mode Execution mode
 * @return boolean if mode matches the indicated modes
 */
bool isSequentialMode(const std::string& mode) {
    std::vector<std::string> sequentialModes = {"sequential", "seq", "s"};
    return std::find(sequentialModes.begin(), sequentialModes.end(), mode) != sequentialModes.end();
}

/**
 * check for parallel mode choices from terminal
 * @param mode Execution mode
 * @return boolean if mode matches the indicated modes
 */
bool isParallelMode(const std::string& mode) {
    std::vector<std::string> sequentialModes = {"parallel", "par", "p"};
    return std::find(sequentialModes.begin(), sequentialModes.end(), mode) != sequentialModes.end();
}

/**
 * Entry point of the program
 */
int main(int argc, char* argv[]) {
    struct timeval t1, t2;
    double sequentialTime = 0.0, parallelTime = 0.0;

    const size_t stSize = strtol(argv[1], NULL, 10);
    string mode = argv[2];

	// optional argument
	int threads = 1; // Default to 1 thread (sequential)

	// Validate command line arguments
	// Check if the number of threads argument is provided and if we are in parallel or both mode
    if ((isParallelMode(mode) || mode == "both") && argc == 4) {
        threads = strtol(argv[3], NULL, 10);
        if (threads < 2) {
            cerr << "Invalid number of threads. Please ensure the number of threads is greater than 2.\n";
            return EXIT_FAILURE;
        }
    }

    if (stSize < 2) {
        cerr << "Invalid array size. Please ensure array size is greater than 2.\n";
        return EXIT_FAILURE;
    }
    
    
    int *data = (int*) aligned_alloc(CACHE_LINE, stSize * sizeof(int));
    int *tmp = (int*) aligned_alloc(CACHE_LINE, stSize * sizeof(int));
    int *ref = (int*) aligned_alloc(CACHE_LINE, stSize * sizeof(int));

    
    if (!data || !tmp || !ref) {
        printf("Memory allocation failed!\n");
        return EXIT_FAILURE;
    }

    printf("Initialization...\n");

    // Initialize array with random values
    srand(95);  
    for (size_t idx = 0; idx < stSize; ++idx) {
        data[idx] = (int)(stSize * (double(rand()) / RAND_MAX));
    }
    copy(data, data + stSize, ref);

    // Calculate and display array size in MiB
    double dSize = (stSize * sizeof(int)) / 1024.0 / 1024.0;
    printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

    /*
        Execution options:
            1. Sequentially: run only original sequential merge sort. by sending run argument (e.g. s).
            2. Parallel: run only added parallel implementation for merge sort. by sending run argument (e.g. p).
            3. run sequential and parallel together and calculate speedup for comparison. by sending run argument "both".
    */

    // Run sequential sort 
    if (isSequentialMode(mode) || mode == "both") {
        gettimeofday(&t1, NULL);
        MsSequential(data, tmp, true, 0, stSize);
        gettimeofday(&t2, NULL);
        sequentialTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        sequentialTime /= 1000.0;
        printf("Sequential execution completed in %f sec.\n", sequentialTime);
    }

    // Run parallel sort 
    if (isParallelMode(mode) || mode == "both") {
        copy(ref, ref + stSize, data);  // Reset array to original state
        gettimeofday(&t1, NULL);
        MsParallel(data, tmp, stSize, threads);
        gettimeofday(&t2, NULL);
        parallelTime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        parallelTime /= 1000.0;
        printf("Parallel execution completed in %f sec.\n", parallelTime);
    }

    // Calculate speedup for both modes
    if (mode == "both") {
        calculateSpeedup(sequentialTime, parallelTime);
    }

    // Verify sorting result
    printf("Verification...");
    if (isSorted(ref, data, stSize)) {
        printf(" successful.\n");
    } else {
        printf(" FAILED.\n");
    }

    free(data);
    free(tmp);
    free(ref);
    return EXIT_SUCCESS;
}