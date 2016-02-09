serial: serial.c
	gcc -o $@ $< -std=c99 -Wall -Wextra -Wpedantic -O3

run: serial
	./serial

.PHONY: clean
clean:
	rm serial
