SHELL = /bin/sh
CC = g++
#WARNINGS = -Wall -Wextra -pedantic -Wshadow
CFLAGS = -funroll-loops -O3 -DNDEBUG -march=native -std=c++17 -pthread
#CFLAGS = -g2 -std=c++17 -pthread

all: bm-text-to-lz77

bm-text-to-lz77:
	$(CC) $(CFLAGS) $(WARNINGS) -o bm-compression ./src/*.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f bm-compression *.o
