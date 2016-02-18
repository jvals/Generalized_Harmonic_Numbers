#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Prototypes
double compute_sum(double v[], uint16_t);

int main(int argc, char *argv[]) {
    // ----- Get input
    if(argc < 2) {
	printf("Usage: ./only_omp.c limit\n limit is an integer\n");
	return 1;
    }

    int nprocs, my_rank;
    uint16_t n;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Rank 0 does I/O
    if(my_rank == 0) {
        uint8_t k = atoi(argv[1]);
	n = pow(2, k);
    }

    // Every rank needs to know the size of the problem
    MPI_Bcast(&n, 1, MPI_UINT16_T, 0, MPI_COMM_WORLD);

    // ----- Generate vector v
    // Vector to hold the partial sums
    double v[n];
    for(uint16_t i = 1; i <= n; i++) {
        v[i] = 1 / (double)(i * i);
    }

    // ----- Compute sum S(n)
    double total_sum = compute_sum(v, n);
    printf("Total Sum: %f\n", total_sum);

    // ----- Compute the error |S - S(n)|
    // The actual sum (computed by wolfram alpha) with 16 digits of accuracy
    double wolfram_sum = 1.644934066848226;
    double error = fabs(wolfram_sum - total_sum);

    // ----- Print out the error |S - S(n)|
    printf("n = %d,\terror = %f\n", n, error);
    
    return 0;
}

double compute_sum(double v[], uint16_t n) {
    double sum = 0.0;
    for(uint16_t i = 1; i <= n; i++) {
        sum += v[i];
    }
    return sum;
}
