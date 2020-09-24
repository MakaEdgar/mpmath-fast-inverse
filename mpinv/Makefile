CC=g++
CFLAGS=-c -fPIC
DIR_OUT=bin
INCLUDE=-I/mnt/b/libs

all: main

main: main.o libmpinv.a
	$(CC) $(DIR_OUT)/main.o -L$(DIR_OUT) -lmpinv -o $(DIR_OUT)/main

main.o: makedir main.cpp
	$(CC) $(CFLAGS) main.cpp -o $(DIR_OUT)/main.o

libmpinv.a: makedir mpinv.cpp
	$(CC) $(INCLUDE) $(CFLAGS) mpinv.cpp -o $(DIR_OUT)/mpinv.o
	ar rcs $(DIR_OUT)/libmpinv.a $(DIR_OUT)/mpinv.o

makedir:
	mkdir -p $(DIR_OUT)

clean:
	rm -rf $(DIR_OUT)/*
