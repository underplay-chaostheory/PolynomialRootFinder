#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complexe.h"

# define PI           3.14159265358979323846

void print_complexe(complexe* z){
    printf("%f +i%f", z->re, z->im);
}

complexe* new_complexe_ref(double a, double b){
    complexe* z = (complexe*)malloc(sizeof(complexe));
    z->re = a;
    z->im = b;
    return z;
}

complexe* new_complexe_polar_ref(double m, double phi){
    complexe* z = (complexe*)malloc(sizeof(complexe));
    z->re = m * cos(phi);
    z->im = m * sin(phi);
    return z;
}

void free_complexe(complexe* z){
    free(z);
}

void free_complexe_array(complexe** array, int array_length){
    for(int i = 0; i < array_length; i++){
        free_complexe(array[i]);
    }
    free(array);
}

complexe* copy_complexe(complexe* z){
    complexe* Z = new_complexe_ref(z->re, z->im);
    return Z;
}
void reset_complexe(complexe* z){
    z->re = 0;
    z->im = 0;
}
void replace_complexe(complexe* z, complexe* r){
    z->re = r->re;
    z->im = r->im;
}

double mod(complexe* z){
    double m = sqrt(z->re*z->re + z->im*z->im);
    return m; 
}

double arg(complexe* z){
    if (z->re == 0){
        if (z->im > 0){
            return PI / 2;
        } else {
            return - PI / 2;
        }
    } else {
        return tan(z->im / z->re);
    }
}

double mod2(complexe* z){
    double m = mod(z);
    return m*m;
}

complexe* conjugate(complexe* z){
    return new_complexe_ref(z->re, -z->im);
}

void conjugate_ip(complexe* z){
    z->im = -z->im;
}

complexe* opposite(complexe* z){
    return new_complexe_ref(-z->re, -z->im);
}

void opposite_ip(complexe* z){
    z->re = -z->re;
    z->im = -z->im;
}

complexe* add(complexe* z1, complexe* z2){
    return new_complexe_ref(z1->re + z2->re, z1->im + z2->im);
}

void add_ip(complexe* z1, complexe* z2){
    z1->re += z2->re;
    z1->im += z2->im;
}

complexe* sub(complexe* z1, complexe* z2){
    return new_complexe_ref(z1->re - z2->re, z1->im - z2->im);
}

void sub_ip(complexe* z1, complexe* z2){
    z1->re -= z2->re;
    z1->im -= z2->im;
}

complexe* mult(complexe* z1, complexe* z2){
    double RE = z1->re*z2->re - z1->im*z2->im;
    double IM = z1->re*z2->im + z1->im*z2->re;
    return new_complexe_ref(RE, IM);
}

void mult_ip(complexe* z1, complexe* z2){
    double RE = z1->re*z2->re - z1->im*z2->im;
    double IM = z1->re*z2->im + z1->im*z2->re;
    z1->re = RE;
    z1->im = IM;
}

void mult_scalar_ip(complexe* z, double x){
    z->re = z->re * x;
    z->im = z->im * x;
    return;
}

complexe* invert(complexe* z){
    complexe* a = conjugate(z);
    double m = mod2(z);
    mult_scalar_ip(z, 1 / m);
    return a;
}

void invert_ip(complexe* z){
    double m = mod2(z);
    conjugate_ip(z);
    z->re = z->re / m;
    z->im = z->im / m;
}

complexe* divide(complexe* z1, complexe* z2){
    double m = mod2(z2);
    complexe* temp = conjugate(z2);
    complexe* a = mult(z1, temp);
    a->re = a->re / m;
    a->im = a->im / m;
    free_complexe(temp);
    return a;
}

void divide_ip(complexe* z1, complexe* z2){
    double m = mod2(z2);
    complexe* temp = conjugate(z2);
    mult_ip(z1, temp);
    z1->re = z1->re / m;
    z1->im = z1->im / m;
    free_complexe(temp);
}

complexe* complexe_square_root(complexe* z){
    double m = mod(z);
    m = sqrt(m);
    double a = arg(z);
    a = a / 2;
    complexe* res = new_complexe_ref(m*cos(a), m*sin(a));
    return res;
}

double distance(complexe* z1, complexe* z2){
    complexe* temp = sub(z1, z2);
    double d = mod(temp);
    free_complexe(temp);
    return d;
}