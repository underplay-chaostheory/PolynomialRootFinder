#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#include "lib/complexe.h"
#include "lib/polynomial.h"

# define PI                 3.1415926535897932384


complexe** complexe_zeros(int len){
    complexe** tab = (complexe**)malloc(len * sizeof(complexe*));
    for(int i = 0; i < len; i++){
        tab[i] = new_complexe_ref(0,0);
    }
    return tab;
}

complexe** set(int degree, int* size, double R, int Njump){
    int nb_circles = (int)ceil(0.26632 * log((double)degree));
    int N = (int)ceil(8.32547 * degree * log((double)degree));

    double* pts_module = (double*)malloc(nb_circles * sizeof(double));
    double K = 1 + sqrt(2);
    double a = (degree-1)/degree;
    
    for(int i = 1; i <= nb_circles; i++){
        pts_module[i-1] = pow(a, (2*i-1)/(4*nb_circles)) * K * R;
    }

    double* pts_argument = (double*)malloc(N * sizeof(double));
    for(int k = 0; k < Njump; k++){
        for(int j = k; j < N; j=j+Njump){
            pts_argument[j] = 2 * PI * j / N;
        }
    }
    
    complexe** Sd = (complexe**)malloc(N * nb_circles * sizeof(complexe*));
    for(int x = 0; x < N; x++){
        for(int y = 0; y < nb_circles; y++){
            Sd[y + x*nb_circles] = new_complexe_polar_ref(pts_module[y], pts_argument[x]);
        }
    }

    free(pts_module);
    free(pts_argument);
    *size = N * nb_circles;
    return Sd;
}

int add_to_set(complexe** roots, complexe* z, int nb_root_found, double EPS){
    int j = 0;
    while (j < nb_root_found){
        if(distance(z, roots[j]) <= EPS){
            break;
        } else {
            j++;
        }
    }
    if(j == nb_root_found){
        replace_complexe(roots[nb_root_found], z);
        return nb_root_found + 1;
    }else{
        return nb_root_found;
    }
}

complexe* pn(int n, complexe* z) {
    //P(n+2) = 2XP(n+1) - P(n)
    complexe* p0 = new_complexe_ref(1.0, 0.0);
    if (n == 0) return p0;

    complexe* p1 = copy_complexe(z);
    if (n == 1) {
        free_complexe(p0);
        return p1;
    }

    complexe *next = NULL;
    for (int i = 2; i <= n; ++i) {
        complexe* temp = mult(z, p1);          
        mult_scalar_ip(temp, 2.0);             
        next = sub(temp, p0);
        free_complexe(temp);

        free_complexe(p0);
        p0 = p1;
        p1 = next;
    }
    free_complexe(p0);
    return next;
}

complexe* dpn(int n, complexe* z){
    //P'(n+2) = 2P(n+1) + 2XP'(n+1) - P'(n)
    if (n == 0){
        complexe* dp0 = new_complexe_ref(0.0, 0.0);
        return dp0;
    }
    if(n == 1){
        complexe* dp1 = new_complexe_ref(1.0, 0.0);
        return dp1;
    }
    if (n == 2){
        complexe* dp2 = copy_complexe(z);
        mult_scalar_ip(z, 4.);
        return dp2;
    }

    complexe** p = complexe_zeros(n);
    complexe** dp = complexe_zeros(n+1);
    complexe* unite = new_complexe_ref(1,0);
    replace_complexe(p[0], unite);
    replace_complexe(dp[1], unite);
    replace_complexe(p[1], z);
    replace_complexe(dp[2], z);
    mult_scalar_ip(dp[2], 4.);
    for(int i = 3; i <= n; i++){
        complexe* temp = mult(p[i-2], z);          
        mult_scalar_ip(temp, 2.0);             
        complexe* next = sub(temp, p[i-3]);
        replace_complexe(p[i-1], next);
        free_complexe(next);
        free_complexe(temp);

        temp = mult(dp[i-1], z);
        mult_scalar_ip(temp, 2.);
        sub_ip(temp, dp[i-2]);
        replace_complexe(dp[i], p[i-1]);
        mult_scalar_ip(dp[i], 2.);
        add_ip(dp[i], temp);
        free_complexe(temp);
    }
    complexe* pz = copy_complexe(dp[n]);
    free_complexe(unite);
    free_complexe_array(p, n);
    free_complexe_array(dp, n+1);
    return pz;
}

double evaluate_mod(int n, complexe* z){
    complexe* pz = pn(n, z);
    double m = mod(pz);
    free_complexe(pz);
    return m;
}

void step(complexe* z0, complexe* prev, int n, double w){
    replace_complexe(prev,z0);

    complexe* z = pn(n, prev);
    complexe* temp = dpn(n, prev);
    divide_ip(z, temp);

    mult_scalar_ip(z, w);

    sub_ip(z, prev);
    opposite_ip(z);

    replace_complexe(z0,z);
    free_complexe(z);
    free_complexe(temp);
}

int root_finder(int n, FILE* OUTPUT_FILE, int ITERATION_LOOP, double EPS, double relaxation){
    time_t t1 = clock();
    int nm_iteration = 0;
    
    int nb_root = n;
    complexe** roots = complexe_zeros(nb_root);
    int nb_root_found = 0;

    //Creating the set Sd
    int Sd_size;
    complexe** Sd = set(n, &Sd_size, 1.0, n);

    complexe** aux = complexe_zeros(Sd_size);
    int completed = 0;
    int set_size = Sd_size;
    int aux_size = 0;
    complexe** set = Sd;
    
    int max_nm_iteration = n * 15000;
    double tol = EPS / n;
    while (nb_root_found < nb_root && completed < Sd_size && nm_iteration < max_nm_iteration){
        int i = 0;
        while(i < set_size && nb_root_found < nb_root && completed < Sd_size){
            complexe* z0 = copy_complexe(set[i]);
            complexe* prev = new_complexe_ref(0,0);

            int k = 0;
            do {
                step(z0, prev, n, relaxation);

                nm_iteration++;
                k++;
            } while (k < ITERATION_LOOP && distance(z0, prev) > tol);

            if(nb_root_found < nb_root && distance(z0, prev) <= tol){
                completed++;
                nb_root_found = add_to_set(roots, z0, nb_root_found, EPS);
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

    time_t t2 = clock();
    double duration = (double)(t2-t1)/1000000.;

    fprintf(OUTPUT_FILE,"%d;%d;%lf\n",n,nb_root_found,duration);
    if(nb_root_found == n){
        for(int k =0; k < n; k++){
            fprintf(OUTPUT_FILE,"%.12f\n", roots[k]->re);
        }
    }
    
    free_complexe_array(roots, nb_root);
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
        int n;
        int k;
        fscanf(INPUT_FILE, "%d\n", &n);
        fscanf(INPUT_FILE, "%d\n", &k);

        if(root_finder(n, OUTPUT_FILE, k, EPS, relaxation) != 0){
            printf("There seems to be an issue.\n");
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