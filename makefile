.PHONY: all
all: main
main: main.c
	gcc -pthread -lm main.c -o main
