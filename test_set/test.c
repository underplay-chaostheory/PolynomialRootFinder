#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "complexe.h"
#include "polynomial.h"
#include "get_test.h"


int root_finder(char* TEST_FILE){
    test* t = get_test_file(TEST_FILE);

    for(int test = 0; test < 1; test ++){
        polynomial* p = get_polynomial(t);
        polynomial* dp = derivate(p);

        polynomial* simple_roots = modulus_polynomials(p, dp);
        print_polynomial(simple_roots);
        free_polynomial(simple_roots);
        simple_roots = pgcd_polynomials(p, dp);
        print_polynomial(simple_roots);
        free_polynomial(simple_roots);
        simple_roots = division_polynomials(p, dp);
        print_polynomial(simple_roots);
        free_polynomial(simple_roots);

        free_polynomial(dp);
        free_polynomial(p);
    }

    close_test(t);
    return EXIT_SUCCESS;
}

int main(int argc, char** argv){
    if(argc != 2){
        return EXIT_FAILURE;
    }
    FILE* INPUT_FILE = fopen(argv[1], "r");
    if(INPUT_FILE == NULL){
        return EXIT_FAILURE;
    }
    printf("Opening %s\n", argv[1]);

    int nb_test;
    char output_file[100];
    double EPS;
    double relaxation;
    fscanf(INPUT_FILE, "%d", &nb_test);
    fscanf(INPUT_FILE, "%s\n", output_file);
    fscanf(INPUT_FILE, "%lf", &EPS);
    fscanf(INPUT_FILE, "%lf", &relaxation);

    printf("%d tests :\n", nb_test);
    for(int i = 0; i < nb_test; i++){
        char test_file[100];
       
        int k;
        fscanf(INPUT_FILE, "%s\n", test_file);
        fscanf(INPUT_FILE, "%d\n", &k);

        if(root_finder(test_file) != 0){
            printf("There seems to be an issue with the test file %s.\n", test_file);
            return EXIT_FAILURE;
        }
    }
    
    if(fclose(INPUT_FILE) != 0){
        printf("There seems to be an issue with input_file.\n");
        return EXIT_FAILURE;
    }
    printf("Closing %s\n", argv[1]);
    return EXIT_SUCCESS;
}