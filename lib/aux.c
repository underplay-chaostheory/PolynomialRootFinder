#include "aux.h"

double radius_complex_hull(polynomial* p){
    //Goershgorin (p must be monic)
    double upper_bound = mod(p->coefficients[0]);
    for(int i = 1; i < p->degree; i++){
        double x = 1 + mod(p->coefficients[i]);
        if(x > upper_bound){
            upper_bound = x;
        }
    }
    return upper_bound + 1;
    //max(1, upperbound) in reality
}

polynomial* simplify(polynomial* p){
    polynomial* dp = derivate(p);
    polynomial* multiple_degree_zeros = pgcd_polynomials(p, dp);
    polynomial* simple_roots = division_polynomials(p, multiple_degree_zeros);
    free_polynomial(dp);
    free_polynomial(multiple_degree_zeros);
    return simple_roots;
}

complexe** complexe_zeros(int len){
    complexe** tab = (complexe**)malloc(len * sizeof(complexe*));
    for(int i = 0; i < len; i++){
        tab[i] = new_complexe_ref(0,0);
    }
    return tab;
}

complexe** set_Njump(int degree, int* size, double R, int Njump){
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

complexe** set_default(int degree, int* size, double R){
    return set_Njump(degree, size, R, 1);
}

complexe** set(int degree, int* size, double R){
    return set_Njump(degree, size, R,1);
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

void step(complexe* z0, complexe* prev, polynomial* p, polynomial* dp, double w){
    replace_complexe(prev,z0);

    complexe* z = evaluate(p, prev);
    complexe* temp = evaluate(dp, prev);
    divide_ip(z, temp);

    mult_scalar_ip(z, w);

    sub_ip(z, prev);
    opposite_ip(z);

    replace_complexe(z0,z);
    free_complexe(z);
    free_complexe(temp);
}

double evaluate_mod(polynomial* p, complexe* z){
    complexe* pz = evaluate(p, z);
    double m = mod(pz);
    free_complexe(pz);
    return m;
}

void solve_deg2(polynomial* p, complexe** roots, int* nb_root){
    if(p->degree == 2){
        complexe* b2 = mult(p->coefficients[1], p->coefficients[1]);
        complexe* temp = copy_complexe(p->coefficients[0]);
        mult_scalar_ip(temp, 4);
        complexe* delta2 = sub(b2, temp);
        free_complexe(b2);
        free_complexe(temp);

        complexe* delta = complexe_square_root(delta2);
        free_complexe(delta2);

        complexe* roots0 = opposite(p->coefficients[1]);
        add_ip(roots0, delta);
        mult_scalar_ip(roots0, 0.5);
        replace_complexe(roots[0], roots0);
        complexe* roots1 = opposite(p->coefficients[1]);
        sub_ip(roots1, delta);
        mult_scalar_ip(roots1, 0.5);
        replace_complexe(roots[1], roots1);

        free_complexe(roots0);
        free_complexe(roots1);
        free_complexe(delta);
        *nb_root = *nb_root + 2;
    } else if (p->degree == 1) {
        complexe* roots0 = opposite(p->coefficients[0]);
        divide_ip(roots0, p->coefficients[1]);
        replace_complexe(roots[0], roots0);
        free_complexe(roots0);
        *nb_root = *nb_root + 1;
    }
}