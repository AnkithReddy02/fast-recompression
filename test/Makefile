SHELL = /bin/sh
CC = g++
#WARNINGS = -Wall -Wextra -pedantic -Wshadow
CFLAGS = -funroll-loops -O3 -DNDEBUG -march=native -std=c++17 -pthread
#CFLAGS = -g2 -std=c++17 -pthread

all: recomp

recomp:
	$(CC) $(CFLAGS) $(WARNINGS) -o test ../src/recompression_definitions.cpp ./src/*.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f recomp *.o
