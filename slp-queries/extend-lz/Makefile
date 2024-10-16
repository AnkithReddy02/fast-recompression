SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x -pthread
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x -pthread

all: extend

extend:
	$(CC) $(CFLAGS) -o extend ./src/main.cpp ./src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f extend *.o
