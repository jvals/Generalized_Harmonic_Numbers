K := 14

serial: serial.c
	gcc -o $@ $< -std=c99 -Wall -Wextra -Wpedantic

run: serial
	./serial $(K)

.PHONY: clean
clean:
	rm serial
