PROCESSES := 4
K := 14

only_mpi: only_mpi.c
	mpicc -o $@ $< -std=c99 -Wall -Wextra -Wpedantic

run: only_mpi
	mpirun -np $(PROCESSES) ./only_mpi $(K)

.PHONY: clean

clean:
	rm only_mpi
