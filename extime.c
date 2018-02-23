/* C code to exactly calculate the expected extinction time for the SIS model on networks

Developed by Petter Holme, September 2017

The input to the program is the network information: First, the number of nodes N, then a list of the edges (assuming the vertices are labeled by numbers from 0 to N - 1). For example, to analyze a triangle, do as follows:

$ ./extime 3 0 1 1 2 2 0

 1 2 4, (2*x^2+3*x+3)/3
 3 5 6, (4*x^2+8*x+9)/6
 7, (4*x^2+8*x+11)/6

The output should be interpreted as follows:

To the left of the comma is a set of automorphically equivalent states. To the right is the expression for the expected extinction time. These will always be a fraction of polynomials---the numerator with a degree N - 1 larger than the degree of the denominator.

*/

#include "poly.h"

GLOBAL g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine eliminates states that are automorphically equivalent
// (so by a relabeling the graph, one can go from one to another)

void get_automorphisms (AUTOCLASS *au) {
	int i, j, me, you;
	igraph_bool_t iso;

	g.autoclass = malloc(g.ns * sizeof(int));

	for (me = 0; me < g.ns; me++) g.autoclass[me] = -1;

	for (me = 1, g.nauto = 0; me < g.ns; me++) if (g.autoclass[me] < 0) {
		g.autoclass[me] = g.nauto++;
		// constructing the color vector corresponding to state me
		for (i = 0; i < g.n; i++)
			VECTOR(g.cme)[i] = (INFECTIOUS(me, i)) ? 1 : 0;
		for (you = me + 1; you < g.ns; you++) if (g.autoclass[you] < 0) {
			// constructing the color vector corresponding to state you
			for (i = 0; i < (int) igraph_vcount(&g.g); i++)
				VECTOR(g.cyou)[i] = (INFECTIOUS(you, i)) ? 1 : 0;
			// checking if the colorings are isomorphic
			igraph_isomorphic_vf2(&g.g, &g.g, &g.cme, &g.cyou, NULL, NULL,
				&iso, NULL, NULL, NULL, NULL, NULL);
			if (iso) g.autoclass[you] = g.autoclass[me];
		}
	}

	// construct the nodes of the state network, one per automorphic equivalence class

	for (i = 0; i < g.nauto; i++) {
		au[i].v = malloc(STATE * sizeof(int));
		au[i].v[N] = 0;
		au[i].p = malloc((g.nauto + 1) * sizeof(fmpz_poly_t));
		for (j = 0; j <= g.nauto; j++) fmpz_poly_init(au[i].p[j]);
	}
	for (me = 1; me < g.ns; me++) au[g.autoclass[me]].v[N]++;
	for (i = 0; i < g.nauto; i++) {
		au[i].v = (int *) realloc(au[i].v, (au[i].v[N] + STATE) * sizeof(int));
		au[i].v[N] = 0;
	}
	for (me = 1; me < g.ns; me++) {
		i = g.autoclass[me];
		au[i].v[STATE + au[i].v[N]++] = me;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void set_up_equations (AUTOCLASS *au) {
	int i, j, ime, iyou, now, *nb, ns2i, ni2s;
	igraph_vs_t vs_me, vs_you;
	igraph_vit_t me, you;

	nb = malloc(g.nauto * sizeof(int));

	igraph_vs_all(&vs_me);

	// for non-absorbing states
	for (i = 0; i < g.nauto; i++) { // if i is the representative member of a class
		// check nodes to generate possible next states
		now = au[i].v[STATE];

		for (j = 0; j < g.nauto; j++) nb[j] = 0;
		ns2i = ni2s = 0;

		igraph_vit_create(&g.g, vs_me, &me);
		
		while (!IGRAPH_VIT_END(me)) {
			ime = (int) IGRAPH_VIT_GET(me);
			if (INFECTIOUS(now, ime)) {
				// positive numbers represent change to S
				// skip transitions to the absorbing state
				j = g.autoclass[RECOVER(now, ime)];
				if (j >= 0) nb[j]++; // ignoring transitions to the absorbing state
				else ni2s++;
				igraph_vs_adj(&vs_you, ime, IGRAPH_ALL);
				igraph_vit_create(&g.g, vs_you, &you);
				while (!IGRAPH_VIT_END(you)) {
					// negative numbers represent change to I
					iyou = (int) IGRAPH_VIT_GET(you);
					if (SUSCEPTIBLE(now, iyou))
						nb[g.autoclass[INFECT(now, iyou)]]--;
					IGRAPH_VIT_NEXT(you);
				}
				igraph_vit_destroy(&you);
			}
			IGRAPH_VIT_NEXT(me);
		}

		igraph_vit_destroy(&me);

		// sum up to normalize probabilities
		for (j = 0; j < g.nauto; j++) {
			if (nb[j] < 0) ns2i -= nb[j];
			else ni2s += nb[j];
		}

		// set the denominator nu
		fmpz_poly_zero(g.a);
		fmpz_poly_set_coeff_si(g.a, 1, -ns2i);
		fmpz_poly_set_coeff_si(g.a, 0, -ni2s);

		// construct the set of equations
		for (j = 0; j < g.nauto; j++) {
			if (nb[j] < 0) { // if S->I
				fmpz_poly_set_coeff_si(au[i].p[j], 1, -nb[j]);
			} else { // if I->S
				if (nb[j] > 0) fmpz_poly_set_coeff_si(au[i].p[j], 0, nb[j]);
				else fmpz_poly_zero(au[i].p[j]); // if no change
			}
		}
		fmpz_poly_set(au[i].p[i], g.a); // the self contribution
		fmpz_poly_set_si(au[i].p[g.nauto], -1); // the constant term
	}

	free(nb);
	free(g.autoclass);

	igraph_vs_destroy(&vs_me);
	igraph_vs_destroy(&vs_you);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// help routine to shift pivot in the gaussian elimination 

void swap_autoclass (int me, int you, AUTOCLASS *a) {
	int *vtmp;
	fmpz_poly_t *ptmp;

	vtmp = a[me].v;
	a[me].v = a[you].v;
	a[you].v = vtmp;

	ptmp = a[me].p;
	a[me].p = a[you].p;
	a[you].p = ptmp;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void solve_equations (AUTOCLASS *au) {
	int i, j, k, jpiv;

	// Gaussian elimination

	for (i = 0; i < g.nauto; i++) {
		for (j = jpiv = i; j < g.nauto; j++)
			if (!fmpz_poly_is_zero(au[i].p[j])) {
				jpiv = j;
				break;
			}

		// swap rows jpiv and i
		if (jpiv > i) swap_autoclass(i, jpiv, au);

		// Gaussian elimination

		for (j = i + 1; j < g.nauto; j++) {
			for (k = i + 1; k <= g.nauto; k++) {
				fmpz_poly_mul(g.a, au[i].p[k], au[j].p[i]);
				fmpz_poly_mul(g.b, au[j].p[k], au[i].p[i]);
				fmpz_poly_sub(au[j].p[k], g.b, g.a);
			}
			fmpz_poly_zero(au[j].p[i]);
			fmpz_poly_clear(au[j].p[i]);

			// simplify
			fmpz_poly_zero(g.a);
			for (k = i + 1; k <= g.nauto; k++)
				if (!fmpz_poly_is_zero(au[j].p[k])) {
					if (fmpz_poly_is_zero(g.a)) fmpz_poly_set(g.a, au[j].p[k]);
					else fmpz_poly_gcd(g.a, g.a, au[j].p[k]);
				}
			for (k = i + 1; k <= g.nauto; k++)
				if (!fmpz_poly_is_zero(au[j].p[k]))
					fmpz_poly_div(au[j].p[k], au[j].p[k], g.a);
		}
	}

	// back substitution

	for (i = g.nauto - 1; i >= 0; i--) {
		for (j = i - 1; j >= 0; j--) {
			// using i,i to delete j,i
			fmpz_poly_zero(g.c);
			for (k = j; k <= g.nauto; k++) if (k != i) {
				fmpz_poly_mul(g.a, au[i].p[k], au[j].p[i]);
				fmpz_poly_mul(g.b, au[j].p[k], au[i].p[i]);
				fmpz_poly_sub(au[j].p[k], g.b, g.a);

				if (!fmpz_poly_is_zero(au[j].p[k])) {
					if (fmpz_poly_is_zero(g.c)) fmpz_poly_set(g.c, au[j].p[k]);
					else fmpz_poly_gcd(g.c, g.c, au[j].p[k]);
				}
			}
			fmpz_poly_clear(au[j].p[i]);

			// simplify
			for (k = j; k <= g.nauto; k++) if (k != i)
				fmpz_poly_div(au[j].p[k], au[j].p[k], g.c);
		}

		// let the denominator have a positive leading coefficient
		if (fmpz_sgn(fmpz_poly_lead(au[i].p[i])) < 0) {
			fmpz_poly_scalar_mul_si(au[i].p[i], au[i].p[i], -1);
			fmpz_poly_scalar_mul_si(au[i].p[g.nauto], au[i].p[g.nauto], -1);
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling i/o

int main (int argc, char *argv[]) {
	int i, j;
	AUTOCLASS *au;

	if (argc < 3) {
		fprintf(stderr, "usage: ./extime [# nodes] [links]\n");
		return 1;
	}

	fmpz_poly_init(g.a);
	fmpz_poly_init(g.b);
	fmpz_poly_init(g.c);

	// reading and setting up the network
	i = atoi(argv[1]);
	g.ns = 1 << i;
	igraph_empty(&g.g, i, IGRAPH_UNDIRECTED);
	igraph_vector_int_init(&g.cme, i);
	igraph_vector_int_init(&g.cyou, i);

	j = argc / 2 - 1;
	for (i = 0; i < j; i++)
		igraph_add_edge(&g.g, atoi(argv[2 + 2 * i]), atoi(argv[3 + 2 * i]));

	g.n = igraph_vcount(&g.g);

	au = malloc(g.ns * sizeof(AUTOCLASS));

	get_automorphisms(au);
	set_up_equations(au);

	igraph_destroy(&g.g);
	igraph_vector_int_destroy(&g.cme);
	igraph_vector_int_destroy(&g.cyou);

	solve_equations(au);
	
	for (i = 0; i < g.nauto; i++) {
		for (j = 0; j < au[i].v[N]; j++) printf(" %d", au[i].v[STATE + j]);
		printf(", ");
		j = fmpz_poly_is_one(au[i].p[i]);
		if (!j) printf("(");
		fmpz_poly_print_pretty(au[i].p[g.nauto], "x");
		if (!j) {
			printf(")/");
			j = (fmpz_poly_degree(au[i].p[i]) > 0);
			if (j) printf("(");
			fmpz_poly_print_pretty(au[i].p[i], "x");
			if (j) printf(")");
		}
		printf("\n");
	}

	for (i = 0; i < g.nauto; i++) {
		fmpz_poly_clear(au[i].p[g.nauto]);
		free(au[i].v);
		free(au[i].p);
	}

	free(au);
	fmpz_poly_clear(g.a);
	fmpz_poly_clear(g.b);
	fmpz_poly_clear(g.c);

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
