#ifndef GET_TEST
#define GET_TEST

#include <stdio.h>
#include <stdlib.h>
#include "complexe.h"
#include "polynomial.h"

struct test{
    FILE* test_file;
    int nb_test;
    int pol_degree;
    int nb_get;
};
typedef struct test test;

test* get_test_file(char* file_name);
void close_test(test* t);
polynomial* get_polynomial(test* t);

#endif