#!/bin/bash

CC = g++ -std=c++11
CFLAGS = -O3
EXEC = search convert
OBJ1 = gdb.o graph_entropy.o

all: $(EXEC)

search: $(OBJ1) main.o
	$(CC) -o $@ $^ $(CFLAGS)

convert: convert.cpp
	$(CC) $(CFLAGS) -o $@ $<
	
%.o: %.cpp %.hpp
	$(CC) -o $@ -c $< $(CFLAGS)
	
%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm $(EXEC) *.o	
