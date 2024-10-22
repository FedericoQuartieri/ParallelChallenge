/*
*  This file is part of Christian's OpenMP software lab 
*
*  Copyright (C) 2016 by Christian Terboven <terboven@itc.rwth-aachen.de>
*  Copyright (C) 2016 by Jonas Hahnfeld <hahnfeld@itc.rwth-aachen.de>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*
*/

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
#include <chrono>
#include <cstdlib> // Per eseguire comandi di sistema come GNUplot



/**
  * helper routine: check if array is sorted correctly
  */
bool isSorted(int ref[], int data[], const size_t size){
	std::sort(ref, ref + size);
	for (size_t idx = 0; idx < size; ++idx){
		if (ref[idx] != data[idx]) {
			return false;
		}
	}
	return true;
}


/**
  * sequential merge step (straight-forward implementation)
  */
// TODO: cut-off could also apply here (extra parameter?)
// TODO: optional: we can also break merge in two halves
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
  * sequential MergeSort
  */
// TODO: remember one additional parameter (depth)
// TODO: recursive calls could be taskyfied
// TODO: task synchronization also is required
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

void MsParallel(int *array, int *tmp, bool inplace, long begin, long end, int cutoff, int depth) {
	if (depth < cutoff) {
		const long half = (begin + end) / 2;


		#pragma omp parallel num_threads(2)
		{
			#pragma omp single
			{
				#pragma omp task shared(array, tmp)
				MsParallel(array, tmp, !inplace, begin, half, cutoff, depth + 1);

				#pragma omp task shared(array, tmp)
				MsParallel(array, tmp, !inplace, half, end, cutoff, depth + 1);

				#pragma omp taskwait
			}
		}

		if (inplace) {
			MsMergeSequential(array, tmp, begin, half, half, end, begin);
		} else {
			MsMergeSequential(tmp, array, begin, half, half, end, begin);
		}
	} else {
		MsSequential(array, tmp, inplace, begin, end);
	}
}


/**
  * Serial MergeSort
  */
// TODO: this function should create the parallel region
// TODO: good point to compute a good depth level (cut-off)

double parallel(const size_t stSize, const size_t cutoff){
	struct timeval t1, t2;
	double etime;

	int *data = (int*) malloc(stSize * sizeof(int));
	int *tmp = (int*) malloc(stSize * sizeof(int));
	int *ref = (int*) malloc(stSize * sizeof(int));

	printf("Initialization...\n");

	srand(95);
	for (size_t idx = 0; idx < stSize; ++idx){
		data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
	}
	std::copy(data, data + stSize, ref);

	double dSize = (stSize * sizeof(int)) / 1024 / 1024;
	printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);


	omp_set_max_active_levels(cutoff+1);
	omp_set_num_threads(pow(2,cutoff));
	//int max_thread = omp_get_max_threads();
	//std::cout << max_thread << std::endl;
	int max_thread = omp_get_max_threads();
	std::cout << "number of threads: " << max_thread << std::endl;
	std::cout << "number of levels: " << omp_get_max_active_levels() << std::endl;
	std::cout << "number of cores: " << omp_get_num_procs() << std::endl;

	gettimeofday(&t1, NULL);
	MsParallel(data, tmp, true, 0, stSize, cutoff, 0);
	//std::sort(data, data + stSize);
	gettimeofday(&t2, NULL);
	etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
	etime = etime / 1000;*/

	printf("Parallel done, took %f sec. Verification...", etime);
	if (isSorted(ref, data, stSize)) {
		printf(" successful.\n");
	}
	else {
		printf(" FAILED.\n");
	}


	free(data);
	free(tmp);
	free(ref);

	return etime;

}


double serial(const size_t stSize){

	struct timeval t1, t2;
	double etime;

	int *data = (int*) malloc(stSize * sizeof(int));
	int *tmp = (int*) malloc(stSize * sizeof(int));
	int *ref = (int*) malloc(stSize * sizeof(int));


	
	for (size_t idx = 0; idx < stSize; ++idx){
		data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
	}
	std::copy(data, data + stSize, ref);


	auto start = std::chrono::high_resolution_clock::now();
		MsSequential(data, tmp, true, 0, stSize);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	/*gettimeofday(&t1, NULL);
	MsSequential(data, tmp, true, 0, stSize);
	gettimeofday(&t2, NULL);
	etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
	etime = etime / 1000;*/

	printf("Sequential done, took %f sec. Verification...", etime);
	if (isSorted(ref, data, stSize)) {
		printf(" successful.\n");
	}
	else {
		printf(" FAILED.\n");
	}


	free(data);
	free(tmp);
	free(ref);
	

	return etime;
}

/** 
  * @brief program entry point
  */
int main(int argc, char* argv[]) {
	// variables to measure the elapsed time


	// expect one command line arguments: array size
	if (argc < 3) {
		printf("Usage: MergeSort.exe <array size> \n");
		printf("\n");
		return EXIT_FAILURE;
	}


	double max = 0;
	int max_cutoff = -1;
	//int cutoff = static_cast<int>(atoi(argv[2]))-1;
	int cutoff = 0;
	for (cutoff; cutoff <  static_cast<int>(atoi(argv[2])); cutoff++) {
		double ptime = parallel(static_cast<int>(atoi(argv[1])), cutoff);
		double stime = serial(static_cast<int>(atoi(argv[1])));
		double speedup = stime/ptime;
		if (speedup > max){
			max = speedup;
			max_cutoff = cutoff;
		}
		std::cout << "speedup with cutoff " << cutoff << " : " << stime/ptime << std::endl<<std::endl;
	}

	std::cout << "max cutoff: " << max_cutoff << "   with speedup: " << max << std::endl;








	std::ofstream file("tempi.csv");
    file << "Array Size,Serial Time,Parallel Time\n";

    // Prova con array di dimensioni variabili
    //std::vector<int> sizes = {1000, 10000, 100000, 500000, 1000000}; // Puoi modificare le dimensioni a tuo piacere
	std::vector<int> sizes;
    for (int i = 1; i <= 10000000; i+=1000000) {
        sizes.push_back(i);
    }
    for (int size : sizes) {
        

        // Misura il tempo per l'algoritmo seriale
        auto start = std::chrono::high_resolution_clock::now();
        serial(size);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time = end - start;

        // Misura il tempo per l'algoritmo parallelo
        start = std::chrono::high_resolution_clock::now();
        parallel(size, static_cast<int>(atoi(argv[2])));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

        // Stampa e salva i risultati
        std::cout << "Array size: " << size << " | Serial: " << serial_time.count() << " s | Parallel: " << parallel_time.count() << " s\n";
        file << size << "," << serial_time.count() << "," << parallel_time.count() << "\n";
    }

    file.close();
	 // Chiama gnuplot per creare il grafico
    std::ofstream gnuplot_file("plot_commands.gp");
    gnuplot_file << "set terminal png size 800,600\n";
    gnuplot_file << "set output 'tempi_grafico.png'\n";
    gnuplot_file << "set title 'Confronto Tempi Seriali vs Paralleli'\n";
    gnuplot_file << "set xlabel 'Dimensione Array'\n";
    gnuplot_file << "set ylabel 'Tempo (s)'\n";
    gnuplot_file << "set grid\n";
	gnuplot_file << "set xrange [0:*]\n";
	gnuplot_file << "set yrange [0:*]\n";
	gnuplot_file << "set datafile separator ','\n";
    gnuplot_file << "plot 'tempi.csv' using 1:2 with lines title 'Seriale', 'tempi.csv' using 1:3 with lines title 'Parallelo'\n";
    gnuplot_file.close();

    // Esegui il comando gnuplot
    system("gnuplot plot_commands.gp");

    std::cout << "Grafico salvato come 'tempi_grafico.png'\n";

    return 0;






	std::ofstream file("tempi.csv");
    file << "Array Size,Serial Time,Parallel Time\n";

    // Prova con array di dimensioni variabili
    //std::vector<int> sizes = {1000, 10000, 100000, 500000, 1000000}; // Puoi modificare le dimensioni a tuo piacere
	std::vector<int> sizes;
    for (int i = 1; i <= 10000000; i+=1000000) {
        sizes.push_back(i);
    }
    for (int size : sizes) {
        

        // Misura il tempo per l'algoritmo seriale
        auto start = std::chrono::high_resolution_clock::now();
        serial(size);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_time = end - start;

        // Misura il tempo per l'algoritmo parallelo
        start = std::chrono::high_resolution_clock::now();
        parallel(size, static_cast<int>(atoi(argv[2])));
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> parallel_time = end - start;

        // Stampa e salva i risultati
        std::cout << "Array size: " << size << " | Serial: " << serial_time.count() << " s | Parallel: " << parallel_time.count() << " s\n";
        file << size << "," << serial_time.count() << "," << parallel_time.count() << "\n";
    }

    file.close();
	 // Chiama gnuplot per creare il grafico
    std::ofstream gnuplot_file("plot_commands.gp");
    gnuplot_file << "set terminal png size 800,600\n";
    gnuplot_file << "set output 'tempi_grafico.png'\n";
    gnuplot_file << "set title 'Confronto Tempi Seriali vs Paralleli'\n";
    gnuplot_file << "set xlabel 'Dimensione Array'\n";
    gnuplot_file << "set ylabel 'Tempo (s)'\n";
    gnuplot_file << "set grid\n";
	gnuplot_file << "set xrange [0:*]\n";
	gnuplot_file << "set yrange [0:*]\n";
	gnuplot_file << "set datafile separator ','\n";
    gnuplot_file << "plot 'tempi.csv' using 1:2 with lines title 'Seriale', 'tempi.csv' using 1:3 with lines title 'Parallelo'\n";
    gnuplot_file.close();

    // Esegui il comando gnuplot
    system("gnuplot plot_commands.gp");

    std::cout << "Grafico salvato come 'tempi_grafico.png'\n";

    return 0;

//-----------------------------------------------------------------------

//---------------------------------




}
