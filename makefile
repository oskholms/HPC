FLAGS= -lm -Ofast -march=native -pthread

.PHONY: all
all: main
main: main.c
	gcc $(FLAGS) main.c -o main
