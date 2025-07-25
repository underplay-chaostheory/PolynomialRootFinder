#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "lib/complexe.h"
#include "lib/polynomial.h"
#include "lib/get_test.h"
#include "lib/aux.h"

int all_found(polynomial* p, polynomial* dp, complexe** roots, int nb_root_found, int nb_root, int ITERATION_LOOP, double EPS, double w, int* nm_iteration_tot){
    //Creating the set Sd
    int Sd_size;
    complexe** Sd = set(p->degree, &Sd_size, 1.0);

    complexe** aux = complexe_zeros(Sd_size);
    int completed = 0;
    int set_size = Sd_size;
    int aux_size = 0;
    complexe** set = Sd;
    
    int nm_iteration = 0;
    int max_nm_iteration = p->degree * 15000;
    double tol = EPS / p->degree;
    while (nb_root_found < nb_root && completed < Sd_size && nm_iteration < max_nm_iteration){
        int i = 0;
        while(i < set_size && nb_root_found < nb_root && completed < Sd_size){
            complexe* z0 = copy_complexe(set[i]);
            complexe* prev = new_complexe_ref(0,0);

            int k = 0;
            do {
                step(z0, prev, p, dp, w);

                nm_iteration++;
                k++;
            } while (k < ITERATION_LOOP && distance(z0, prev) > tol);
            
            if(nb_root_found < nb_root && distance(z0, prev) > tol){
                completed++;
                nb_root_found = add_to_set(roots, z0, nb_root_found, 2*EPS);
            } else {
                replace_complexe(aux[aux_size], z0);
                aux_size ++;
            }
            free_complexe(z0);
            free_complexe(prev);
            i++;
        }
        set_size = aux_size;
        aux_size = 0;
        complexe** temp = set;
        set = aux;
        aux = temp;
    }

    free_complexe_array(set, Sd_size);
    free_complexe_array(aux, Sd_size);
    *nm_iteration_tot = *nm_iteration_tot + nm_iteration;
    return nb_root_found;
}

int root_finder(char* TEST_FILE, FILE* OUTPUT_FILE, int ITERATION_LOOP, double EPS, double relaxation){
    //Relaxed Newton's Method for w =/= 1
    double w = relaxation;

    //Open test file
    test* t = get_test_file(TEST_FILE);

    for(int test = 0; test < t->nb_test; test ++){
        polynomial* p = get_polynomial_v2(t);

        if(p->degree < 3){
            printf("Erreur dans get_test : %d\n", p->degree);
            free_polynomial(p);
            continue;
        }

        time_t t1 = clock();
        int nm_iteration = 0;

        //Note that P and Q = P / pgdc(P,P') have same roots (without multiplicity)
        polynomial* simple_roots = simplify(p);
        polynomial* dsr = derivate(simple_roots);
        if(simple_roots->degree < 3){
            printf("Erreur dans simplify\n");
            free_polynomial(p);
            free_polynomial(simple_roots);
            free_polynomial(dsr);
            continue;
        }
        
        //We apply the algorithm; if not all roots found, we divide P by (X - r1)*..*(X - rk)
        //where r1, .., rk are the roots that have been found
        int nb_root = simple_roots->degree;
        complexe** roots = complexe_zeros(nb_root);
        int nb_root_found = 0;

        while (nb_root_found < nb_root && simple_roots->degree >= 3){
            int new_nb_root_found = all_found(simple_roots, dsr, roots, nb_root_found, nb_root, ITERATION_LOOP, EPS, w, &nm_iteration);
            
            if(new_nb_root_found == nb_root){
                nb_root_found = new_nb_root_found;
                break;
            }

            int found = new_nb_root_found - nb_root_found;
            if(found == 0){
                break;
            }
            polynomial* root_found_prod = developp(&roots[nb_root_found], found);
            polynomial* r = division_polynomials(simple_roots, root_found_prod);
            
            free_polynomial(root_found_prod);
            free_polynomial(simple_roots);
            free_polynomial(dsr);

            simple_roots = r;
            dsr = derivate(simple_roots);
            nb_root_found = new_nb_root_found;
        }
        if(simple_roots->degree < 3 && nb_root_found != nb_root){
            solve_deg2(simple_roots, &roots[nb_root_found], &nb_root_found);
        }

        //Freeing the variables
        free_polynomial(simple_roots);
        free_polynomial(dsr);

        time_t t2 = clock();
        double duration = (double)(t2-t1)/1000000.;

        if (nb_root_found == nb_root){
            fprintf(OUTPUT_FILE,"%d;%d;%d;%lf\n",ITERATION_LOOP,nb_root,nm_iteration,duration);
            for(int i = 0; i < nb_root; i++){
                fprintf(OUTPUT_FILE, "%.15lf+i%.12lf\n", roots[i]->re, roots[i]->im);
            }
        } else {
            fprintf(OUTPUT_FILE, "-1\n");
        }
        
        free_complexe_array(roots, nb_root);
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

    FILE* OUTPUT_FILE = fopen(output_file, "w");
    if(OUTPUT_FILE == NULL){
        return EXIT_FAILURE;
    }
    printf("Opening %s\n", output_file);

    printf("%d tests :\n", nb_test);
    for(int i = 0; i < nb_test; i++){
        char test_file[100];
       
        int k;
        fscanf(INPUT_FILE, "%s\n", test_file);
        fscanf(INPUT_FILE, "%d\n", &k);

        fprintf(OUTPUT_FILE, "New test file\n");
        fprintf(OUTPUT_FILE, "%s\n", test_file);

        if(root_finder(test_file, OUTPUT_FILE, k, EPS, relaxation) != 0){
            printf("There seems to be an issue with the test file %s.\n", test_file);
            return EXIT_FAILURE;
        }
    }
    
    if(fclose(INPUT_FILE) != 0){
        printf("There seems to be an issue with input_file.\n");
        return EXIT_FAILURE;
    }
    printf("Closing %s\n", argv[1]);
    if(fclose(OUTPUT_FILE) != 0){
        printf("There seems to be an issue with output_file.\n");
        return EXIT_FAILURE;
    }
    printf("Closing %s\n", output_file);
    return EXIT_SUCCESS;
}
