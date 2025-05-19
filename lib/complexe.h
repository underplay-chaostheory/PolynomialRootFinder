#ifndef COMPLEXE_H
#define COMPLEXE_H

struct complexe {
  double re;
  double im;
};
typedef struct complexe complexe;

void print_complexe(complexe* z);

complexe* new_complexe_ref(double a, double b);
complexe* new_complexe_polar_ref(double m, double phi);
void free_complexe(complexe* z);
void free_complexe_array(complexe** array, int array_length);
complexe* copy_complexe(complexe* z);
void reset_complexe(complexe* z);
void replace_complexe(complexe* z, complexe* r);

double mod(complexe* z);
double arg(complexe* z);
double mod2(complexe* z);

complexe* conjugate(complexe* z);
void conjugate_ip(complexe* z);
complexe* opposite(complexe* z);
void opposite_ip(complexe* z);
complexe* add(complexe* z1, complexe* z2);
void add_ip(complexe* z1, complexe* z2);
complexe* sub(complexe* z1, complexe* z2);
void sub_ip(complexe* z1, complexe* z2);
complexe* mult(complexe* z1, complexe* z2);
void mult_ip(complexe* z1, complexe* z2);
void mult_scalar_ip(complexe* z1, double x);
complexe* invert(complexe* z);
void invert_ip(complexe* z);
complexe* divide(complexe* z1, complexe* z2);
void divide_ip(complexe* z1, complexe* z2);
complexe* complexe_square_root(complexe* z);

double distance(complexe* z1, complexe* z2);

#endif