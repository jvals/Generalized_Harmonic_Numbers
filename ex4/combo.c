#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

// Prototypes
double compute_sum(double v[], uint16_t elements_per_process);

int main(int argc, char *argv[]) {
    // Check input
    if (argc < 2) {
        printf("Usage: ./only_omp.c limit\n limit is an integer\n");
        return 1;
    }

    int nprocs, my_rank, provided;
    uint16_t n;

    // Initialize MPI. MPI_THREAD_FUNNELED specifies that only the main thread
    // will make MPI calls. This allows MPI to do some optimization.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // We need to check if MPI is providing the level of thread support
    // that we require.
    if (provided != MPI_THREAD_FUNNELED) {
        printf("MPI_THREAD_FUNNELED is not available.\n");
        return 1;
    }

    // Rank 0 does I/O
    if (my_rank == 0) {
        uint8_t k = atoi(argv[1]);
        n = pow(2, k);
    }

    // Every rank needs to know the size of the problem
    MPI_Bcast(&n, 1, MPI_UINT16_T, 0, MPI_COMM_WORLD);

    // Generate vector v
    // Vector to hold the partial sums
    double *v;
    if (my_rank == 0) {
        v = (double *) malloc(n * sizeof(double));
        for (uint16_t i = 1; i <= n; i++) {
            v[i-1] = 1 / (double)(i * i);
        }
    }

    // Each process only needs to allocate memory for it's own chunk of data
    uint16_t elements_per_process = (n/nprocs);
    double *scattered_v = (double *) malloc(elements_per_process * sizeof(double));
    
    // Send parts of the generated vector to every rank
    MPI_Scatter(v,elements_per_process,MPI_DOUBLE, 
        scattered_v, elements_per_process,MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compute sum S(n)
    double total_sum;
    double partial_sum = compute_sum(scattered_v, elements_per_process);
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

    if (my_rank == 0) {
        free(v);
    }
    free(scattered_v);
    MPI_Finalize();
    return 0;
}

double compute_sum(double v[], uint16_t elements_per_process) {
    double sum = 0.0;
#   pragma omp parallel for reduction(+:sum) num_threads(4) schedule(guided)
    for (uint16_t i = 0; i < elements_per_process; i++) {
        sum += v[i];
    }

    return sum;
}
