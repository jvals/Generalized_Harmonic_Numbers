omp: only_omp.c
	/usr/local/bin/gcc-5 -o $@ $< -Wall -Wextra -Wpedantic -fopenmp

run: omp
	./omp

.PHONY: clean
clean:
	rm omp
