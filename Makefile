SRC_DIR=src
BIN_DIR=bin

CPP=g++
CPPFLAGS=-c -fPIC -DNDEBUG -fno-math-errno -O3 -march=native -std=c++11
INCLUDE=-Ilibs

all: main

main: main.o libmpinv.a
	$(CPP) $(BIN_DIR)/main.o -L$(BIN_DIR) -lmpinv -o $(BIN_DIR)/main

main.o: makedir $(SRC_DIR)/main.cpp
	$(CPP) $(INCLUDE) $(CPPFLAGS) $(SRC_DIR)/main.cpp -o $(BIN_DIR)/main.o

libmpinv.a: makedir $(SRC_DIR)/mpinv.cpp
	$(CPP) $(INCLUDE) $(CPPFLAGS) $(SRC_DIR)/mpinv.cpp -o $(BIN_DIR)/mpinv.o
	ar rcs $(BIN_DIR)/libmpinv.a $(BIN_DIR)/mpinv.o

makedir:
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)/*
