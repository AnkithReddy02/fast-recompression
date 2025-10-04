SHELL = /bin/sh
CC = g++
#WARNINGS = -Wall -Wextra -pedantic -Wshadow
CFLAGS = -funroll-loops -O3 -DNDEBUG -march=native -std=c++17 -pthread
#CFLAGS = -g2 -std=c++17 -pthread

INCLUDES = -Iinclude

all: recomp

recomp:
	$(CC) $(CFLAGS) $(INCLUDES) $(WARNINGS) -o recomp ./src/*.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f recomp *.o
