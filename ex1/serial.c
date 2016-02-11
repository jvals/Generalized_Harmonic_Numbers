#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#define MAX 16384 // 2^14

// Prototypes
double compute_sum(double v[], uint16_t);

int main() {
    // ----- Generate vector v
    // Vector to hold the partial sums
    double v[MAX];
    for(uint16_t i = 1; i <= MAX; i++) {
        v[i] = 1 / (double)(i * i);
    }

    // ----- Compute sum S(n)
    double total_sum = compute_sum(v, MAX);
    printf("Total Sum: %f\n", total_sum);

    // ----- Compute the error |S - S(n)|
    // The actual sum (computed by wolfram alpha) with 16 digits of accuracy
    double wolfram_sum = 1.644934066848226;
    // New vector containing the errors
    double e[11];
    for(uint8_t k = 3; k <= 14; k++) {
        double partial_sum = compute_sum(v, pow(2,k));
        e[k-3] = fabs(wolfram_sum - partial_sum);
    }

    // ----- Print out the error |S - S(n)|
    for(uint8_t i = 0; i <= 11; i++) {
        printf("k = %d,\terror = %f\n", i+3, e[i]);
    }

    return 0;
}

double compute_sum(double v[], uint16_t n) {
    double sum = 0.0;
    for(uint16_t i = 1; i <= n; i++) {
        sum += v[i];
    }
    return sum;
}
