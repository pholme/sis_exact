SRC = .
FLINT_DIR = "/usr/local/include/flint"
IGRAPH_DIR = "/usr/local/include/igraph"
CFLAGS = -W -Wall -I$(FLINT_DIR) -I$(IGRAPH_DIR) -Ofast -march=native -fomit-frame-pointer
LDFLAGS = -lflint -lmpfr -lgmp -lpthread -ligraph
CC = gcc

OBJ1 = o/extime.o

all : extime

extime : $(OBJ1)
	$(CC) -o $@ $^ $(LDFLAGS)

o/extime.o : $(SRC)/extime.c $(SRC)/poly.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/extime.c -o $@
