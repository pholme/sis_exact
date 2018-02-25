// header files for the code to calculate SIS extinction times by Petter Holme

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fmpz_poly.h>
#include <igraph.h>

#define N 0
#define NINF 1  // this is not needed for this version of the program, but is an example of how one can store info about an equivalence class
#define SMALL 2 // .. same for this ..
#define LARGE 3 // and this (so you can delete these three and set STATE to 1)
#define STATE 4

#define INFECTIOUS(x,y) ((x) & (1 << (y))) // if node y is I in state x
#define SUSCEPTIBLE(x,y) (!((x) & (1 << (y)))) // if node y is S in state x
#define INFECT(x,y) ((x) | (1 << (y))) // the state x with y infected 
#define RECOVER(x,y) ((x) - (1 << (y))) // the state x with y susceptible

typedef struct { // automorphic equivalence classes of states
	int *v; // vector containing info about autoclass
	fmpz_poly_t *p; // extinction time of states in the class
} AUTOCLASS;

typedef struct {
	int n; // number of nodes
	int ns; // number of SIS states
	int nauto; // number of automorphically distinct states
	int *autoclass;
	fmpz_poly_t a, b, c; //scratch space
	igraph_t g; // the graph
	igraph_vector_int_t cme, cyou; // color vectors
} GLOBAL;
