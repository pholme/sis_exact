This code calculates the expected time to extinction in the SIS model from all possible configurations of susceptible and infectious and any network. It uses graph isomorphism to reduce the computational time if there are symmetries in the underlying network. The code requires the FLINT library http://www.flintlib.org/ for polynomial algebra and igraph http://igraph.org/c/ for the graph isomorphism routine.

The code is compiled by make (it assumes a directory o/ where it can put the object files). You'd probably have to edit the Makefile to reflect your environment.

You can run the program like:

`./extime 3 0 1 1 2`

Where the first argument, 3, is the number of nodes. The rest is a list of links to define the graph, the node indices are assumed to range from 0 to _N_ - 1.

The output is:

```
1 4, (4*x^4+16*x^3+35*x^2+34*x+12)/(16*x^2+28*x+12)
2, (2*x^3+7*x^2+14*x+6)/(8*x+6)
3 6, (4*x^3+16*x^2+35*x+18)/(16*x+12)
5, (4*x^4+20*x^3+53*x^2+52*x+18)/(16*x^2+28*x+12)
7, (4*x^4+20*x^3+57*x^2+62*x+22)/(16*x^2+28*x+12)
```

Before the comma, there is a list of equivalent configurations. After the comma is the expression for the _x_-values of the corresponding configurations.

For more info, see the accompanying paper: Holme & Tupikina, Epidemic extinction in networks: Insights from the 12,110 smallest graphs.
