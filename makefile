.PHONY: all
all: main
main: main.c
	gcc -Ofast -pthread -lm main.c -o main
