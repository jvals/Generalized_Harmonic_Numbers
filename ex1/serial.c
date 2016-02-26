#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

// Prototypes
double compute_sum(double v[], uint16_t);

int main(int argc, char *argv[]) {
    // Check input
    if(argc < 2) {
        printf("Usage: ./serial.c limit\n limit is an integer\n");
        return 1;
    }
    uint8_t k = atoi(argv[1]);
    uint16_t n = pow(2, k);

    // Generate vector v
    // Vector to hold the partial sums
    double v[n];
    for(uint16_t i = 1; i <= n; i++) {
        v[i-1] = 1. / (i * i);
    }

    // Compute sum S(n)
    double total_sum = compute_sum(v, n);
    printf("Total Sum: %.16f\n", total_sum);

    // Compute the error |S - S(n)|
    // The actual sum (computed by wolfram alpha) with 16 digits of accuracy
    double wolfram_sum = 1.644934066848226;
    double error = fabs(wolfram_sum - total_sum);

    // Print out the error |S - S(n)|
    printf("n = %d,\terror = %.16f\n", n, error);
    
    return 0;
}

double compute_sum(double v[], uint16_t n) {
    double sum = 0.0;
    for(uint16_t i = 0; i < n; i++) {
        sum += v[i];
    }
    return sum;
}
