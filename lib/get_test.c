#include "get_test.h"

test* get_test_file(char* file_name){
    FILE* f = fopen(file_name, "r");
    test* t = malloc(sizeof(test));
    if(fscanf(f, "%d\n", &(t->nb_test))){
        t->pol_degree = -1;
        t->nb_get = 0;
        t->test_file = f;
        return t;
    }else{
        t->pol_degree = -1;
        t->nb_get = 0;
        t->test_file = NULL;
        t->nb_test = 0;
        return t;
    }
}

void close_test(test* t){
    if(fclose(t->test_file) == 0){
        free(t);
    } else {
        free(t);
        printf("WARNING !");
    }
}

polynomial* get_polynomial(test* t){
    if(t->nb_get > t->nb_test){
        polynomial* p = new_polynomial(-1);
        t->pol_degree = -1;
        return p;
    }
    if(fscanf(t->test_file, "%d\n", &(t->pol_degree))){
        polynomial* p = new_polynomial(t->pol_degree);
        for(int i = 0; i < t->pol_degree; i++){
            double x;
            double y;
            fscanf(t->test_file, "%lf", &x);
            fscanf(t->test_file, "%lf", &y);
            complexe* z = new_complexe_ref(x,y);
            set_p_i(p, z, i);
            free_complexe(z);
        }
        complexe* unite = new_complexe_ref(1,0);
        set_p_i(p, unite, p->degree);
        free_complexe(unite);

        t->nb_get ++;
        return p;
    }else{
        t->pol_degree = -1;
        return new_polynomial(-1);
    }
}

polynomial* get_polynomial_v2(test* t){
    if(t->nb_get > t->nb_test){
        polynomial* p = new_polynomial(-1);
        t->pol_degree = -1;
        return p;
    }
    if(fscanf(t->test_file, "%d\n", &(t->pol_degree))){
        polynomial* p = new_polynomial(t->pol_degree);
        for(int i = 0; i < t->pol_degree; i++){
            char s[100];
            fgets(s, 100, t->test_file);
        }
        for(int i = 0; i < t->pol_degree; i++){
            double x;
            double y;
            fscanf(t->test_file, "%lf", &x);
            fscanf(t->test_file, "%lf", &y);
            complexe* z = new_complexe_ref(x,y);
            set_p_i(p, z, i);
            free_complexe(z);
        }
        complexe* unite = new_complexe_ref(1,0);
        set_p_i(p, unite, p->degree);
        free_complexe(unite);

        t->nb_get ++;
        return p;
    }else{
        t->pol_degree = -1;
        return new_polynomial(-1);
    }
}
