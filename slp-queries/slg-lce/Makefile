SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x -pthread
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x -pthread

all: test

test:
	$(CC) $(CFLAGS) -o test ./src/main.cpp ./src/utils.cpp ./src/karp_rabin_hashing.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
