#ifndef LIB
#define LIB

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "complexe.h"
#include "polynomial.h"
#include "get_test.h"

# define PI                 3.1415926535897932384


double radius_complex_hull(polynomial* p);
polynomial* simplify(polynomial* p);

complexe** complexe_zeros(int len);
complexe** set(int degree, int* size, double R);
int add_to_set(complexe** roots, complexe* z, int nb_root_found, double EPS);

void step(complexe* z0, complexe* prev, polynomial* p, polynomial* dp, double w);
double evaluate_mod(polynomial* p, complexe* z);

void solve_deg2(polynomial* p, complexe** roots, int* nb_root);

#endif