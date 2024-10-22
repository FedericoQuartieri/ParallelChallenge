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
#include <cstdlib> // Per eseguire comandi di sistema come GNUplot
#include <sstream>

#include <chrono>
#include <thread>



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


//parallel recursive 
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
		//if over cutoff use the sequential algorithm
		MsSequential(array, tmp, inplace, begin, end);
	}
}


//-------------------------------------------------------

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

	std::cout << "number of threads: " << omp_get_max_threads() << std::endl;
	std::cout << "number of levels: " << omp_get_max_active_levels() << std::endl;
	std::cout << "number of cores: " << omp_get_num_procs() << std::endl;

	gettimeofday(&t1, NULL);
	MsParallel(data, tmp, true, 0, stSize, cutoff, 0);
	//std::sort(data, data + stSize);
	gettimeofday(&t2, NULL);
	etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
	etime = etime / 1000;

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

	printf("Initialization...\n");

	srand(95);
	for (size_t idx = 0; idx < stSize; ++idx){
		data[idx] = (int) (stSize * (double(rand()) / RAND_MAX));
	}
	std::copy(data, data + stSize, ref);

	double dSize = (stSize * sizeof(int)) / 1024 / 1024;
	printf("Sorting %zu elements of type int (%f MiB)...\n", stSize, dSize);

	gettimeofday(&t1, NULL);
	MsSequential(data, tmp, true, 0, stSize);
	gettimeofday(&t2, NULL);
	etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
	etime = etime / 1000;

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

	if (argc < 3) {
		printf("Usage: MergeSort.exe <array size> \n");
		printf("\n");
		return EXIT_FAILURE;
	}

	int max_possible_size = atoi(argv[1]);
	int max_possible_cutoff = atoi(argv[2]);
	int h = atoi(argv[3]);
	int divider_for_avg = atoi(argv[4]);


    // Prova con array di dimensioni variabili
    //std::vector<int> sizes = {1000, 10000, 100000, 500000, 1000000}; // Puoi modificare le dimensioni a tuo piacere
	std::vector<int> sizes;
    for (int i = max_possible_size/h; i <= max_possible_size; i+=max_possible_size/h) {
		sizes.push_back(i);
    }

	double max_speedup = 0;
	int max_cutoff = -1;
	int max_size = 0;

	for (int cutoff = 0; cutoff < max_possible_cutoff; cutoff++){
		std::stringstream ss;
		ss << "tempi_size_ " << cutoff << ".csv";
		std::string s = ss.str();

		std::ofstream file(s);
		file << "Array Size,Serial Time,Parallel Time\n";

		for (int size : sizes) {
			double serial_time=0.0;
			double parallel_time=0.0;
			for (int i = 0; i < divider_for_avg; i++) {
 				serial_time += serial(size);
 				parallel_time += parallel(size, cutoff);
			}
			serial_time /= divider_for_avg;
			parallel_time /= divider_for_avg;
			std::cout << "serial time: " << serial_time << std::endl;
			std::cout << "parallel time: " << parallel_time << std::endl;

			double speedup = serial_time/parallel_time;
			
			std::cout << "speedup with cutoff " << cutoff << " : " << speedup << std::endl<<std::endl;
			if (speedup > max_speedup){
				max_speedup = speedup;
				max_cutoff = cutoff;
				max_size = size;
			}

			// Stampa e salva i risultati
			std::cout << "Array size: " << size << " | Serial: " << serial_time << " s | Parallel: " << parallel_time << " s\n\n\n";
			file << size << "," << serial_time << "," << parallel_time << "\n";
		}

		file.close();

		std::stringstream ss1;
		ss1 << "tempi_size_" << cutoff << ".png";
		// Chiama gnuplot per creare il grafico
		std::ofstream gnuplot_file("plot_commands.gp");
		gnuplot_file << "set terminal png size 800,600\n";
		gnuplot_file << "set output '" << ss1.str() << "'\n";//+
		gnuplot_file << "set title 'Confronto Tempi Seriali vs Paralleli'\n";
		gnuplot_file << "set xlabel 'Dimensione Array'\n";
		gnuplot_file << "set ylabel 'Tempo (s)'\n";
		gnuplot_file << "set grid\n";
		gnuplot_file << "set xrange [" << max_possible_size/h << ":*]\n";//+
		gnuplot_file << "set yrange [0:*]\n";
		gnuplot_file << "set datafile separator ','\n";
		gnuplot_file << "plot '"+s+"' using 1:2 with lines title 'Seriale', '"+s+"' using 1:3 with lines title 'Parallelo'\n";
		gnuplot_file.close();

		// Esegui il comando gnuplot
		system("gnuplot plot_commands.gp");

		std::cout << "Grafico salvato\n\n";
		std::cout << "---------------------------------------------\nNew cutoff\n---------------------------------------------\n\n";
	}

	std::ofstream file("tempo_cutoff.csv");
	file << "Cutoff,Serial Time,Parallel Time\n";

    for (int cutoff = 0; cutoff < max_possible_cutoff; cutoff++) {
		
		double serial_time=0.0;
		double parallel_time=0.0;
		for (int i = 0; i < divider_for_avg; i++) {
			serial_time += serial(max_possible_size);
			parallel_time += parallel(max_possible_size, cutoff);
		}
		serial_time /= divider_for_avg;
		parallel_time /= divider_for_avg;
		std::cout << "serial time: " << serial_time << std::endl;
		std::cout << "parallel time: " << parallel_time << std::endl;

		double speedup = serial_time/parallel_time;
		
		std::cout << "speedup with cutoff " << cutoff << " : " << speedup << std::endl<<std::endl;
		if (speedup > max_speedup){
			max_speedup = speedup;
			max_cutoff = cutoff;
		}

		// Stampa e salva i risultati
		std::cout << "Array size: " << max_possible_size << " | Serial: " << serial_time << " s | Parallel: " << parallel_time << " s\n\n\n";
		file << cutoff << "," << serial_time << "," << parallel_time << "\n";
	}

	file.close();


	std::ofstream gnuplot_file("plot_commands.gp");
	gnuplot_file << "set terminal png size 800,600\n";
	gnuplot_file << "set output 'tempi_cutoff_max_size.png'\n";//+
	gnuplot_file << "set title 'Confronto Tempi Seriali vs Paralleli'\n";
	gnuplot_file << "set xlabel 'Cutoff'\n";
	gnuplot_file << "set ylabel 'Tempo (s)'\n";
	gnuplot_file << "set grid\n";
	gnuplot_file << "set xrange [0:*]\n";
	gnuplot_file << "set yrange [0:*]\n";
	gnuplot_file << "set datafile separator ','\n";
	gnuplot_file << "plot 'tempo_cutoff.csv' using 1:2 with lines title 'Seriale', 'tempo_cutoff.csv' using 1:3 with lines title 'Parallelo'\n";
	gnuplot_file.close();

	double dSize = (max_size * sizeof(int)) / 1024 / 1024;
	std::cout << "max cutoff: " << max_cutoff << " with speedup: " << max_speedup << " and " << "size: " << dSize << " MiB (" << max_size << " elements)" << std::endl;


	// Esegui il comando gnuplot
	system("gnuplot plot_commands.gp");




    return 0;
//-----------------------------------------------------------------------

//---------------------------------




}
