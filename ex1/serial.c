#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#define MAX 16384 // 2^14

int main() {
    // Vector to hold the partial sums
    double v[MAX];

    // Sum variable
    double sum = 0;

    // The actual sum (computed by wolfram alpha) with 16 digits of accuracy
    double actual_sum = 1.644934066848226;

    // Set the loop variable here, because we only need to initialize it once
    uint16_t n = 1;
    
    // n = 2 ^ k, k = {3, 4, ... , 14}
    for(uint8_t k = 3; k <= 14; k++) {
	for(; n <= pow(2,k); n++) {
	    // Compute the result and put it in the vector, even though it's
	    // strictly not necessary on a single processor 
	    v[n] = 1 / (double)(n * n);
	    sum += v[n];
    	}
	// Print the iteration, and the error
	printf("k = %d,\t error = %f\n", k, actual_sum - sum);
    }
    
    return 0;
}
