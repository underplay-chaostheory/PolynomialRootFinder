#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "complexe.h"

struct polynomials
{
    int degree;
    complexe** coefficients;
};
typedef struct polynomials polynomial;

polynomial* new_polynomial(int degree);
void free_polynomial(polynomial* p);
polynomial* copy_polynomial(polynomial* p);
void replace_polynomial(polynomial* a, polynomial* b);
void adjust_degree(polynomial* p);

void print_polynomial(polynomial* p);

complexe* evaluate(polynomial* p, complexe* z);

polynomial* derivate(polynomial* p);

void mult_polynomial_by_complexe(polynomial* p, complexe* z);
void mult_polynomial_by_scalar(polynomial* p, double x);
void normalize(polynomial* p);
polynomial* add_polynomials(polynomial* p1, polynomial* p2);
polynomial* sub_polynomials(polynomial* p1, polynomial* p2);
polynomial* mult_polynomials(polynomial* p1, polynomial* p2);
polynomial* modulus_polynomials(polynomial* a, polynomial* b);
polynomial* pgcd_polynomials(polynomial* a, polynomial* b);
polynomial* division_polynomials(polynomial* a, polynomial* b);

polynomial* develop(complexe** roots, int nb_roots);

#endif