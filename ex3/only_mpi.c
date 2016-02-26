#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Prototypes
double compute_sum(double v[], uint16_t n, int my_rank, int nprocs);

int main(int argc, char *argv[]) {
    // Check input
    if (argc < 2) {
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
    if (my_rank == 0) {
        uint8_t k = atoi(argv[1]);
        n = pow(2, k);
    }

    // Every rank needs to know the size of the problem
    MPI_Bcast(&n, 1, MPI_UINT16_T, 0, MPI_COMM_WORLD);

    // Generate vector v
    // Vector to hold the partial sums
    double *v = (double *) malloc(n * sizeof(double));
    if (my_rank == 0) {
        for (uint16_t i = 1; i <= n; i++) {
            v[i-1] = 1 / (double)(i * i);
        }
    }
    // Broadcast the generated vector to every rank
    MPI_Bcast(v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compute sum S(n)
    double total_sum;
    double partial_sum = compute_sum(v, n, my_rank, nprocs);
    MPI_Reduce(&partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("Total Sum: %.16f\n", total_sum);

        // Compute the error |S - S(n)|
        // The actual sum (computed by wolfram alpha) with 16 digits of accuracy
        double wolfram_sum = 1.644934066848226;
        double error = fabs(wolfram_sum - total_sum);

        // Print out the error |S - S(n)|
        printf("n = %d,\terror = %.16f\n", n, error);
    }

    free(v);
    MPI_Finalize();
    return 0;
}

double compute_sum(double v[], uint16_t n, int my_rank, int nprocs) {
    double sum = 0.0;

    uint16_t my_start = (n * my_rank) / nprocs;
    uint16_t my_stop = (n * (my_rank + 1)) / nprocs;

    for (uint16_t i = my_start; i < my_stop; i++) {
        sum += v[i];
    }

    return sum;
}
