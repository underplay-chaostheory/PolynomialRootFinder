#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polynomial.h"
#include "complexe.h"

double LIM = 0.0000000000001;


polynomial* new_polynomial(int degree){
    polynomial* p = (polynomial*)malloc(sizeof(polynomial));
    p->degree = degree;
    if(degree == -1){
        p->coefficients = NULL;
        return p;
    }
    p->coefficients = (complexe**)malloc((degree + 1) * sizeof(complexe*));
    for(int i = 0; i <= p->degree; i++){
        p->coefficients[i] = new_complexe_ref(0,0);
    }
    return p;
}

void free_polynomial(polynomial* p){
    if(p->degree == -1){
        free(p);
        return;
    }
    free_complexe_array(p->coefficients, p->degree + 1);
    free(p);
    return;
}
void polynomial_reset(polynomial* p){
    if(p->degree > -1){
        free_complexe_array(p->coefficients, p->degree + 1);
    }
    p->degree = -1;
    p->coefficients = NULL;
}

polynomial* copy_polynomial(polynomial* p){
    polynomial* copy = new_polynomial(p->degree);
    if(p->degree > -1){
        for(int i = 0; i <= p->degree; i++){
            free_complexe(copy->coefficients[i]);
            copy->coefficients[i] = copy_complexe(p->coefficients[i]);
        }
    }
    return copy;
}

void replace_polynomial(polynomial* a, polynomial* b){
    polynomial_reset(a);
    if(b->degree > -1){
        a->degree = b->degree;
        a->coefficients = (complexe**)malloc((b->degree + 1)* sizeof(complexe*));
        for(int i = 0; i <= a->degree; i++){
            a->coefficients[i] = copy_complexe(b->coefficients[i]);
        }
    }
}

void adjust_degree(polynomial* p){
    while (p->degree >= 0 && mod(p->coefficients[p->degree]) < LIM){
        free_complexe(p->coefficients[p->degree]);
        p->degree --;
    }
    if(p->degree == -1){
        free(p->coefficients);
        p->coefficients = NULL;
    }
}

void print_polynomial(polynomial* p){
    if(p->degree == -1){
        printf("P = 0");
        return;
    }
    for(int k = 0; k <= p->degree; k++){
        if (k == 0){
            print_complexe(p->coefficients[k]);
        } else if (k == 1){
            printf(" + ");
            print_complexe(p->coefficients[k]);
            printf("x");
        } else {
            printf(" + ");
            print_complexe(p->coefficients[k]);
            printf("x^%d", k);
        }
    }
    printf("\n");
    return;
}


complexe* evaluate(polynomial* p, complexe* z){
    if(p->degree == -1){
        complexe* evaluation = new_complexe_ref(0,0);
        return evaluation;
    }
    complexe* evaluation = copy_complexe(p->coefficients[p->degree]);
    for(int i = p->degree - 1; i >= 0; i--){
        mult_ip(evaluation, z);
        add_ip(evaluation, p->coefficients[i]);
    }
    return evaluation;
}

polynomial* derivate(polynomial* p){
    if(p->degree < 1){
        polynomial* dp = new_polynomial(-1);
        return dp;
    }
    polynomial* dp = new_polynomial(p->degree - 1);
    for(int i = 0; i < p->degree; i++){
        complexe* z = copy_complexe(p->coefficients[i+1]);
        z->re *= (i+1);
        z->im *= (i+1);
        free_complexe(dp->coefficients[i]);
        dp->coefficients[i] = z;
    }
    return dp;
}


void mult_polynomial_by_complexe(polynomial* p, complexe* z){
    if(mod(z) < LIM){
        polynomial_reset(p);
        return;
    }
    for(int i = 0; i <= p->degree; i++){
        mult_ip(p->coefficients[i], z);
    }
    return;
}

void mult_polynomial_by_scalar(polynomial* p, double x){
    if(-LIM < x && x < LIM){
        polynomial_reset(p);
        return;
    }
    for(int i = 0; i <= p->degree; i++){
        mult_scalar_ip(p->coefficients[i], x);
    }
    return;
}

void normalize(polynomial* p){
    mult_polynomial_by_scalar(p, 1 / mod(p->coefficients[p->degree]));
    complexe* z = new_complexe_ref(1,0);
    replace_complexe(p->coefficients[p->degree],z);
}

polynomial* add_polynomials(polynomial* p1, polynomial* p2){
    if(p1->degree == -1){
        return copy_polynomial(p2);
    }
    if(p2->degree == -1){
        return copy_polynomial(p1);
    }
    polynomial* small;
    polynomial* big;
    if (p1->degree < p2->degree){
        small = p1;
        big = p2;
    } else {
        small = p2;
        big = p1;
    }
    polynomial* s = new_polynomial(big->degree);
    for(int i = 0; i <= small->degree; i++){
        complexe* z = add(small->coefficients[i], big->coefficients[i]);
        free_complexe(s->coefficients[i]);
        s->coefficients[i] = z;
    }
    for(int i = small->degree + 1; i <= big->degree; i++){
        free_complexe(s->coefficients[i]);
        s->coefficients[i] = copy_complexe(big->coefficients[i]);
    }
    adjust_degree(s);
    return s;
}

polynomial* sub_polynomials(polynomial* p1, polynomial* p2){
    if(p1->degree == -1){
        polynomial* s = copy_polynomial(p2);
        for(int i = 0; i <= s->degree; i++){
            opposite_ip(s->coefficients[i]);
        }
        return s;
    }
    if(p2->degree == -1){
        return copy_polynomial(p1);
    }
    if (p1->degree < p2->degree){
        polynomial* s = new_polynomial(p2->degree);
        for(int i = 0; i <= p1->degree; i++){
            complexe* z = sub(p1->coefficients[i], p2->coefficients[i]);
            free_complexe(s->coefficients[i]);
            s->coefficients[i] = z;
        }
        for(int i = p1->degree + 1; i <= p2->degree; i++){
            free_complexe(s->coefficients[i]);
            s->coefficients[i] = opposite(p2->coefficients[i]);
        }
        adjust_degree(s);
        return s;
    } else {
        polynomial* s = new_polynomial(p1->degree);
        for(int i = 0; i <= p2->degree; i++){
            complexe* z = sub(p1->coefficients[i], p2->coefficients[i]);
            free_complexe(s->coefficients[i]);
            s->coefficients[i] = z;
        }
        for(int i = p2->degree + 1; i <= p1->degree; i++){
            free_complexe(s->coefficients[i]);
            s->coefficients[i] = copy_complexe(p1->coefficients[i]);
        }
        adjust_degree(s);
        return s;
    }
}


polynomial* mult_polynomials(polynomial* p1, polynomial* p2){
    polynomial* p = new_polynomial(p1->degree + p2->degree);
    for(int k = 0; k <= p1->degree; k++){
        for(int l = 0; l <= p2->degree; l++){
            complexe* z = mult(p1->coefficients[k], p2->coefficients[l]);
            add_ip(p->coefficients[k+l], z);
            free_complexe(z);
        }
    }
    return p;
}

polynomial* modulus_polynomials(polynomial* a, polynomial* b) {
    if (b->degree == -1) {
        polynomial* p = new_polynomial(-1);
        return p;
    }

    polynomial* p = copy_polynomial(a);

    while (p->degree >= b->degree) {  
        complexe* z = divide(p->coefficients[p->degree], b->coefficients[b->degree]);

        for (int i = 0; i <= b->degree; i++) {
            complexe* t = mult(b->coefficients[b->degree - i], z);
            sub_ip(p->coefficients[p->degree - i], t);
            free_complexe(t);
        }
        
        free_complexe(z);
        adjust_degree(p);
    }
    return p;
}

polynomial* pgcd_polynomials(polynomial* a, polynomial* b) {
    polynomial* a_copy = copy_polynomial(a);
    polynomial* b_copy = copy_polynomial(b);

    while (b_copy->degree != -1) {
        polynomial* r = modulus_polynomials(a_copy, b_copy);

        replace_polynomial(a_copy, b_copy);
        replace_polynomial(b_copy, r);
        free_polynomial(r);
    }

    free_polynomial(b_copy);
    return a_copy;
}

polynomial* division_polynomials(polynomial* a, polynomial* b) {
    if (b->degree == -1) {
        return new_polynomial(-1);  // Cas où b est nul, division impossible
    }

    if (b->degree == 0) {
        if (mod(b->coefficients[0]) < LIM && mod(b->coefficients[0]) > (1 / LIM)) {
            return new_polynomial(-1);
        } else {
            complexe* z = invert(b->coefficients[0]);
            if (!z) return NULL;

            polynomial* q = copy_polynomial(a);
            if (!q) {
                free_complexe(z);
                return NULL;
            }

            mult_polynomial_by_complexe(q, z);
            free_complexe(z);
            return q;
        }
    }

    if (a->degree < b->degree) {
        return new_polynomial(-1);  // Degré trop bas pour être divisé
    }

    polynomial* quotient = new_polynomial(a->degree - b->degree);  // Stocke le quotient
    if (!quotient) return NULL;

    polynomial* r = copy_polynomial(a);  // Copie du dividende
    if (!r) {
        free_polynomial(quotient);
        return NULL;
    }

    while (r->degree >= b->degree) {
        int deg_diff = r->degree - b->degree;
        complexe* dom = divide(r->coefficients[r->degree], b->coefficients[b->degree]);

        // Stocker le coefficient dans le quotient
        replace_complexe(quotient->coefficients[deg_diff], dom);

        // Soustraction du multiple de b au reste
        for (int i = 0; i <= b->degree; i++) {
            complexe* t = mult(b->coefficients[b->degree - i], dom);
            complexe* u = sub(r->coefficients[r->degree - i], t);
            replace_complexe(r->coefficients[r->degree - i], u);

            free_complexe(t);
            free_complexe(u);
        }

        free_complexe(dom);
        adjust_degree(r);
    }

    free_polynomial(r);
    complexe* unite = new_complexe_ref(1.,0.);
    replace_complexe(quotient->coefficients[quotient->degree], unite);
    free_complexe(unite);
    return quotient;
}

polynomial* mult_polynomial_by_X_mod(polynomial*p, int n){
    polynomial* xp = new_polynomial(n);
    int k;
    if(p->degree < n){
        k = p->degree;
    } else {
        k = n-1;
    }
    for(int i = k; i >= 0; i--){
        replace_complexe(xp->coefficients[i+1], p->coefficients[i]);
    }
    return xp;
}

polynomial* develop(complexe** roots, int nb_roots){
    polynomial* p = new_polynomial(nb_roots);
    polynomial* temp;

    free_complexe(p->coefficients[0]);
    p->coefficients[0] = new_complexe_ref(1, 0);

    for(int k = 0; k < nb_roots; k++){
        temp = mult_polynomial_by_X_mod(p, nb_roots);
        mult_polynomial_by_complexe(p, roots[k]);
        polynomial* q = sub_polynomials(temp, p);
        free_polynomial(temp);
        free_polynomial(p);
        p = q;
    }
    return p;
}