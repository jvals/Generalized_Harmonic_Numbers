PROCESSES := 4
K := 14

combo: combo.c
	mpicc -o $@ $< -std=c99 -fopenmp -Wall -Wextra -Wpedantic

run: combo
	mpirun -np $(PROCESSES) ./combo $(K)

.PHONY: clean

clean:
	rm combo
