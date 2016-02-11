serial: serial.c
	gcc -o $@ $< -std=c99 -Wall -Wextra -Wpedantic

run: serial
	./serial

.PHONY: clean
clean:
	rm serial
