#include <stdio.h>      
#include <stdlib.h>     
#include <errno.h>      
#include <sys/time.h>   

#include <iostream>     
#include <algorithm>    
#include <fstream>      
#include <cstdlib>      
#include <cstdio>       

#include <cmath>        
#include <ctime>        
#include <cstring>      

#include <omp.h>        
#include <cstdlib>      
#include <sstream>      

#include <chrono>       
#include <thread>       

#include <vector>
#include <string>
#include <iomanip> // For std::setprecision   
#include <map>
#include <set>

/**
 * @brief Performs bilinear interpolation on a 2D grid.
 *
 * @param x_values Vector of x-axis values.
 * @param y_values Vector of y-axis values.
 * @param f_values 2D grid of function values.
 * @param x Interpolation point along the x-axis.
 * @param y Interpolation point along the y-axis.
 * @return Interpolated value at point (x, y).
 */
double bilinearInterpolation(
    const std::vector<int>& x_values,
    const std::vector<int>& y_values,
    const std::vector<std::vector<double>>& f_values,
    double x, double y) {

    // Find the indices surrounding the x and y values
    auto it_x = std::lower_bound(x_values.begin(), x_values.end(), x);
    auto it_y = std::lower_bound(y_values.begin(), y_values.end(), y);

    if (it_x == x_values.end() || it_y == y_values.end() || it_x == x_values.begin() || it_y == y_values.begin()) {
        // Out of bounds, return 0 or extrapolate as needed
        return 0.0;
    }

    int x1_idx = std::distance(x_values.begin(), it_x) - 1;
    int x2_idx = x1_idx + 1;
    int y1_idx = std::distance(y_values.begin(), it_y) - 1;
    int y2_idx = y1_idx + 1;

    double x1 = x_values[x1_idx];
    double x2 = x_values[x2_idx];
    double y1 = y_values[y1_idx];
    double y2 = y_values[y2_idx];

    double fQ11 = f_values[x1_idx][y1_idx];
    double fQ21 = f_values[x2_idx][y1_idx];
    double fQ12 = f_values[x1_idx][y2_idx];
    double fQ22 = f_values[x2_idx][y2_idx];

    double denom = (x2 - x1) * (y2 - y1);

    double fxy = (fQ11 * (x2 - x) * (y2 - y) +
                  fQ21 * (x - x1) * (y2 - y) +
                  fQ12 * (x2 - x) * (y - y1) +
                  fQ22 * (x - x1) * (y - y1)) / denom;

    return fxy;
}

bool isSorted(int ref[], int data[], const size_t size){
    // Sort the reference array using std::sort
    std::sort(ref, ref + size);

    // Compare each element of the reference array with the data array
    for (size_t idx = 0; idx < size; ++idx){
        if (ref[idx] != data[idx]) {
            // If any element is different, return false
            return false;
        }
    }

    // All elements are the same, return true
    return true;
}

void MsMergeSequential(int *out, int *in, long begin1, long end1, long begin2, long end2, long outBegin) {
    long left = begin1;
    long right = begin2;
    long idx = outBegin; // Index in the output array

    while (left < end1 && right < end2) {
        if (in[left] <= in[right]) {
            // If the current element in the left subarray is less than or equal to the one in the right subarray
            out[idx] = in[left]; // Copy it to the output array
            left++; // Move to the next element in the left subarray
        } else {
            // If the current element in the right subarray is less
            out[idx] = in[right]; // Copy it to the output array
            right++; // Move to the next element in the right subarray
        }
        idx++; // Move to the next position in the output array
    }

    // Copy any remaining elements from the left subarray
    while (left < end1) {
        out[idx] = in[left];
        left++;
        idx++;
    }

    // Copy any remaining elements from the right subarray
    while (right < end2) {
        out[idx] = in[right];
        right++;
        idx++;
    }
}

void MsSequential(int *array, int *tmp, bool inplace, long begin, long end) {
    if (begin < (end - 1)) {
        // If the subarray has more than one element
        const long half = (begin + end) / 2; // Find the midpoint

        // Recursively sort the first half
        MsSequential(array, tmp, !inplace, begin, half);

        // Recursively sort the second half
        MsSequential(array, tmp, !inplace, half, end);

        // Merge the two sorted halves
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }
    } else if (!inplace) {
        tmp[begin] = array[begin];
    }
}

void MsParallel(int *array, int *tmp, bool inplace, long begin, long end, int cutoff, int depth) {
    if (depth < cutoff) {
        // If the current depth is less than the cutoff, proceed in parallel
        const long half = (begin + end) / 2; 

        #pragma omp parallel num_threads(2) // Start a parallel region with 2 threads
        {
            #pragma omp single
            {
                // Create a task for the left half
                #pragma omp task shared(array, tmp)
                MsParallel(array, tmp, !inplace, begin, half, cutoff, depth + 1);

                // Create a task for the right half
                #pragma omp task shared(array, tmp)
                MsParallel(array, tmp, !inplace, half, end, cutoff, depth + 1);

                // Wait for both tasks to complete
                #pragma omp taskwait
            }
        }

        // Merge the two sorted halves
        if (inplace) {
            MsMergeSequential(array, tmp, begin, half, half, end, begin);
        } else {
            MsMergeSequential(tmp, array, begin, half, half, end, begin);
        }
    } else {
        // If the current depth is equal to or greater than the cutoff, proceed sequentially
        MsSequential(array, tmp, inplace, begin, end);
    }
}

double parallel(const size_t stSize, const size_t cutoff){
    struct timeval t1, t2; // Variables to store start and end times
    double etime; // Variable to store elapsed time
    
    // Allocate memory for the data arrays
    int *data = (int*) malloc(stSize * sizeof(int)); // Array to be sorted
    int *tmp = (int*) malloc(stSize * sizeof(int)); // Temporary array for merging
    int *ref = (int*) malloc(stSize * sizeof(int)); // Reference array for verification

    printf("Initialization...\n");

    // Initialize the data array with random values
    srand(95);
    for (size_t idx = 0; idx < stSize; ++idx){
        data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
    }
    // Copy the data to the reference array for later verification
    std::copy(data, data + stSize, ref);

    // Calculate the size of the data in MiB
    double dSize = (stSize * sizeof(int)) / 1024 / 1024;
    printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

    // Set OpenMP parameters
    omp_set_max_active_levels(cutoff+1); // Set the maximum number of nested parallel regions
    omp_set_num_threads(pow(2,cutoff)); // Set the number of threads based on the cutoff

    // Print OpenMP information
    std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
    std::cout << "Number of levels: " << omp_get_max_active_levels() << std::endl;
    std::cout << "Number of processors: " << omp_get_num_procs() << std::endl;

    // Get the start time
    gettimeofday(&t1, NULL);

    // Perform the parallel MergeSort
    MsParallel(data, tmp, true, 0, stSize, cutoff, 0);

    // Get the end time
    gettimeofday(&t2, NULL);

    // Calculate the elapsed time in seconds
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;

    printf("Parallel done, took %f sec. Verification...", etime);

    // Verify if the data array is sorted correctly
    if (isSorted(ref, data, stSize)) {
        printf(" successful.\n");
    }
    else {
        printf(" FAILED.\n");
    }

    // Free the allocated memory
    free(data);
    free(tmp);
    free(ref);

    // Return the elapsed time
    return etime;
}

double serial(const size_t stSize){
    struct timeval t1, t2; // Variables to store start and end times
    double etime; // Variable to store elapsed time
    
    // Allocate memory for the data arrays
    int *data = (int*) malloc(stSize * sizeof(int)); // Array to be sorted
    int *tmp = (int*) malloc(stSize * sizeof(int)); // Temporary array for merging
    int *ref = (int*) malloc(stSize * sizeof(int)); // Reference array for verification

    printf("Initialization...\n");

    // Initialize the data array with random values
    srand(95); // Seed the random number generator
    for (size_t idx = 0; idx < stSize; ++idx){
        data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
    }
    // Copy the data to the reference array for later verification
    std::copy(data, data + stSize, ref);

    // Calculate the size of the data in MiB
    double dSize = (stSize * sizeof(int)) / 1024 / 1024;
    printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

    // Get the start time
    gettimeofday(&t1, NULL);

    // Perform the sequential MergeSort
    MsSequential(data, tmp, true, 0, stSize);

    // Get the end time
    gettimeofday(&t2, NULL);

    // Calculate the elapsed time in seconds
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;

    printf("Sequential done, took %f sec. Verification...", etime);

    // Verify if the data array is sorted correctly
    if (isSorted(ref, data, stSize)) {
        printf(" successful.\n");
    }
    else {
        printf(" FAILED.\n");
    }

    // Free the allocated memory
    free(data);
    free(tmp);
    free(ref);

    // Return the elapsed time
    return etime;
}

/**
 * @brief Performs interpolation of timing data and writes results to CSV files.
 */
void performInterpolation() {
    // Read existing data from 'tempi_3d.csv' for interpolation
    std::ifstream data_file("tempi_3d.csv");
    if (!data_file.is_open()) {
        std::cerr << "Error opening tempi_3d.csv for reading." << std::endl;
        return;
    }

    std::string line;
    std::getline(data_file, line); // Skip header line

    // Create a map to store times based on cutoff and size
    std::map<int, std::map<int, std::pair<double, double>>> data_map;

    while (std::getline(data_file, line)) {
        std::stringstream ss(line);
        std::string value;
        int cutoff_val, size_val;
        double serial_time_val, parallel_time_val;

        std::getline(ss, value, ',');
        cutoff_val = std::stoi(value);

        std::getline(ss, value, ',');
        size_val = std::stoi(value);

        std::getline(ss, value, ',');
        serial_time_val = std::stod(value);

        std::getline(ss, value, ',');
        parallel_time_val = std::stod(value);

        data_map[cutoff_val][size_val] = {serial_time_val, parallel_time_val};
    }
    data_file.close();

    // Extract unique cutoff values and array sizes
    std::set<int> cutoff_values_set;
    std::set<int> size_values_set;

    for (const auto& cutoff_pair : data_map) {
        cutoff_values_set.insert(cutoff_pair.first);
        for (const auto& size_pair : cutoff_pair.second) {
            size_values_set.insert(size_pair.first);
        }
    }

    std::vector<int> cutoff_values(cutoff_values_set.begin(), cutoff_values_set.end());
    std::vector<int> size_values(size_values_set.begin(), size_values_set.end());

    // Create grids for serial and parallel times
    std::vector<std::vector<double>> parallel_time_grid(cutoff_values.size(), std::vector<double>(size_values.size()));
    std::vector<std::vector<double>> serial_time_grid(cutoff_values.size(), std::vector<double>(size_values.size()));

    for (size_t i = 0; i < cutoff_values.size(); ++i) {
        for (size_t j = 0; j < size_values.size(); ++j) {
            int cutoff = cutoff_values[i];
            int size = size_values[j];
            parallel_time_grid[i][j] = data_map[cutoff][size].second;
            serial_time_grid[i][j] = data_map[cutoff][size].first;
        }
    }

    // Prepare to write the interpolated data to CSV files
    std::ofstream parallel_interpolated_file("tempi_3d_interpolated_parallel.csv");
    std::ofstream serial_interpolated_file("tempi_3d_interpolated_serial.csv");
    parallel_interpolated_file << "Cutoff,Array Size,Parallel Time\n";
    serial_interpolated_file << "Cutoff,Array Size,Serial Time\n";

    // Define the number of points for interpolation
    int num_cutoff_points = 100;
    int num_size_points = 100;

    // Generate interpolated cutoff and size values
    double cutoff_min = cutoff_values.front();
    double cutoff_max = cutoff_values.back();
    double size_min = size_values.front();
    double size_max = size_values.back();

    double cutoff_step = (cutoff_max - cutoff_min) / (num_cutoff_points - 1);
    double size_step = (size_max - size_min) / (num_size_points - 1);

    for (int i = 0; i < num_cutoff_points; ++i) {
        double interp_cutoff = cutoff_min + i * cutoff_step;
        for (int j = 0; j < num_size_points; ++j) {
            double interp_size = size_min + j * size_step;

            // Perform bilinear interpolation for both parallel and serial times
            double interp_parallel_time = bilinearInterpolation(
                cutoff_values, size_values, parallel_time_grid, interp_cutoff, interp_size);

            double interp_serial_time = bilinearInterpolation(
                cutoff_values, size_values, serial_time_grid, interp_cutoff, interp_size);

            // Write the interpolated data to the respective files
            parallel_interpolated_file << std::fixed << std::setprecision(6)
                              << interp_cutoff << ","
                              << interp_size << ","
                              << interp_parallel_time << "\n";

            serial_interpolated_file << std::fixed << std::setprecision(6)
                              << interp_cutoff << ","
                              << interp_size << ","
                              << interp_serial_time << "\n";
        }
    }
    parallel_interpolated_file.close();
    serial_interpolated_file.close();
}

void createGraphSize(int cutoff, int max_possible_size, int h, std::string s){
// Create a gnuplot script to generate a plot for this cutoff
    std::stringstream ss1;
    ss1 << "tempi_size_" << cutoff << ".png"; 
    std::ofstream gnuplot_file("plot_commands.gp");
    gnuplot_file << "set terminal png size 800,600\n";
    gnuplot_file << "set output '" << ss1.str() << "'\n";
    gnuplot_file << "set title 'Comparison of Serial vs Parallel Times with cutoff "<< cutoff << std::endl;
    gnuplot_file << "set xlabel 'Array Size'\n";
    gnuplot_file << "set ylabel 'Time (s)'\n";
    gnuplot_file << "set grid\n";
    gnuplot_file << "set xrange [" << max_possible_size/h << ":*]\n";
    gnuplot_file << "set yrange [0:*]\n";
    gnuplot_file << "set datafile separator ','\n";
    gnuplot_file << "plot '" << s << "' using 1:2 with lines title 'Serial', '" << s << "' using 1:3 with lines title 'Parallel'\n";
    gnuplot_file.close();

    // Execute the gnuplot command to generate the plot
    system("gnuplot plot_commands.gp");

    std::cout << "Plot saved\n\n";
}

void createGraphCutoff(){
    // Create a gnuplot script to generate a plot of time vs cutoff
    std::ofstream gnuplot_file("plot_commands.gp");
    gnuplot_file << "set terminal png size 800,600\n";
    gnuplot_file << "set output 'tempi_cutoff_max_size.png'\n";
    gnuplot_file << "set title 'Comparison of Serial vs Parallel Times'\n";
    gnuplot_file << "set xlabel 'Cutoff'\n";
    gnuplot_file << "set ylabel 'Time (s)'\n";
    gnuplot_file << "set grid\n";
    gnuplot_file << "set xrange [0:*]\n";
    gnuplot_file << "set yrange [0:*]\n";
    gnuplot_file << "set datafile separator ','\n";
    gnuplot_file << "plot 'tempo_cutoff.csv' using 1:2 with lines title 'Serial', 'tempo_cutoff.csv' using 1:3 with lines title 'Parallel'\n";
    gnuplot_file.close();
}


void executeCutoffSizes(int max_possible_size, int max_possible_cutoff, int h, int divider_for_avg){
    
    double max_speedup = 0; // Variable to store the maximum speedup achieved
    int max_cutoff = -1; // Variable to store the cutoff that gives maximum speedup
    int max_size = 0; // Variable to store the array size that gives maximum speedup

    // Create a vector to store different array sizes to test
    std::vector<int> sizes;
    for (int i = max_possible_size/h; i <= max_possible_size; i += max_possible_size/h) {
        sizes.push_back(i); // Add sizes from max_size/h up to max_size, incrementing by max_size/h
    }

    // Open a CSV file to store 3D data (cutoff, array size, times)
    std::ofstream file_3d("tempi_3d.csv");
    file_3d << "Cutoff,Array Size,Serial Time,Parallel Time\n";

    // Loop over different cutoff values
    for (int cutoff = 0; cutoff < max_possible_cutoff; cutoff++){
        // Create a filename for the CSV file corresponding to this cutoff
        std::stringstream ss;
        ss << "tempi_size_" << cutoff << ".csv";
        std::string s = ss.str();

        // Open a CSV file to store the times for different array sizes at this cutoff
        std::ofstream file(s);
        file << "Array Size,Serial Time,Parallel Time\n";

        // Loop over different array sizes
        for (int size : sizes) {
            double serial_time = 0.0; // Variable to accumulate serial times
            double parallel_time = 0.0; // Variable to accumulate parallel times

            // Run the sort multiple times to average the time
            for (int i = 0; i < divider_for_avg; i++) {
                serial_time += serial(size); // Perform serial sort and accumulate time
                parallel_time += parallel(size, cutoff); // Perform parallel sort and accumulate time
            }
            // Calculate the average times
            serial_time /= divider_for_avg;
            parallel_time /= divider_for_avg;

            std::cout << "Serial time: " << serial_time << " s" << std::endl;
            std::cout << "Parallel time: " << parallel_time << " s" << std::endl;

            // Calculate the speedup
            double speedup = serial_time / parallel_time;
            
            std::cout << "Speedup with cutoff " << cutoff << " : " << speedup << std::endl << std::endl;

            // Check if this is the maximum speedup achieved so far
            if (speedup > max_speedup){
                max_speedup = speedup;
                max_cutoff = cutoff;
                max_size = size;
            }

            // Write the data to the 3D CSV file
            file_3d << cutoff << "," << size << "," << serial_time << "," << parallel_time << "\n";

            // Print and save the results
            std::cout << "Array size: " << size << " | Serial: " << serial_time << " s | Parallel: " << parallel_time << " s\n\n\n";
            file << size << "," << serial_time << "," << parallel_time << "\n";
        }

        file.close(); // Close the CSV file for this cutoff


        createGraphSize(cutoff, max_possible_size, h, s);
        
        std::cout << "---------------------------------------------\nNew cutoff\n---------------------------------------------\n\n";
    }

    file_3d.close();

    
    double dSize = (max_size * sizeof(int)) / 1024 / 1024;
    std::cout << "Max cutoff: " << max_cutoff << " with speedup: " << max_speedup << " and size: " << dSize << " MiB (" << max_size << " elements)" << std::endl;


    // Call the interpolation function
    performInterpolation();
    system("python3 graph.py");

}


void executeCutoff(int max_possible_size, int max_possible_cutoff, int h, int divider_for_avg){

    // Now test the effect of different cutoff values on a fixed maximum array size
    std::ofstream file("tempo_cutoff.csv");
    file << "Cutoff,Serial Time,Parallel Time\n";
    
    for (int cutoff = 0; cutoff < max_possible_cutoff; cutoff++) {
        double serial_time = 0.0; // Variable to accumulate serial times
        double parallel_time = 0.0; // Variable to accumulate parallel times

        // Run the sort multiple times to average the time
        for (int i = 0; i < divider_for_avg; i++) {
            serial_time += serial(max_possible_size);
            parallel_time += parallel(max_possible_size, cutoff);
        }
        // Calculate the average times
        serial_time /= divider_for_avg;
        parallel_time /= divider_for_avg;

        std::cout << "Serial time: " << serial_time << " s" << std::endl;
        std::cout << "Parallel time: " << parallel_time << " s" << std::endl;

        // Calculate the speedup
        double speedup = serial_time / parallel_time;
        
        std::cout << "Speedup with cutoff " << cutoff << " : " << speedup << std::endl << std::endl;

        
        // Print and save the results
        std::cout << "Array size: " << max_possible_size << " | Serial: " << serial_time << " s | Parallel: " << parallel_time << " s\n\n\n";
        file << cutoff << "," << serial_time << "," << parallel_time << "\n";
    }

    file.close();

    createGraphCutoff();
    // Execute the gnuplot command to generate the plot
    system("gnuplot plot_commands.gp");
}


/** 
 * @brief Program entry point.
 */
int main(int argc, char* argv[]) {

    if (argc < 5) {
        printf("Usage: MergeSort.exe <array size> <max cutoff> <h> <divider_for_avg>\n");
        printf("\n");
        return EXIT_FAILURE;
    }

    // Parse command-line arguments
    int max_possible_size = atoi(argv[1]); // Maximum array size
    int max_possible_cutoff = atoi(argv[2]); // Maximum cutoff depth
    int h = atoi(argv[3]); // Used to calculate the array size increments
    int divider_for_avg = atoi(argv[4]); // Number of runs to average the time

    executeCutoff(max_possible_size, max_possible_cutoff, h, divider_for_avg);
    executeCutoffSizes(max_possible_size, max_possible_cutoff, h, divider_for_avg);

    return 0;
}
