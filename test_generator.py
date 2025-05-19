from random import*
from sys import*
from math import*

#file format :
#     number of polynomials
#     degree of P0
#     Re(P0_1)
#     IM(PO_1)
#     Re(P0_2)
#     .
#     .
#     .
#     degree of P1
#     Re(P1_1)
#     ...

#     The degree-th coef ins't shown, as it is equal to 1

class Complex:
    def __init__(self, a, b):
        self.re = a*cos(b)
        self.im = a*sin(b)
    
    def mult_scalar(self, r):
        return Complex(self.re*r, self.im*r)

    def oppose(self):
        self.re = - self.re
        self.im = - self.im
    
    def add(self, z):
        if (type(z) != Complex):
            ValueError('Add parameter must be a Complex')
        else:
            return Complex(self.re + z.re, self.im + z.im)
    
    def sub(self, z):
        if (type(z) != Complex):
            ValueError('Sub parameter must be a Complex')
        else:
            return Complex(self.re - z.re, self.im - z.im)
        
    def mult(self, z):
        if (type(z) != Complex):
            ValueError('Mult parameter must be a Complex')
        else:
            return Complex(self.re*z.re - self.im*z.im, self.im*z.re + self.re*z.im)
    
    def mult_ip(self, z):
        re = self.re*z.re - self.im*z.im
        im = self.im*z.re + self.re*z.im
        self.re = re
        self.im = im

    def mod(self):
        m = self.re ** 2 + self.im ** 2
        return m
    
    def copy(self):
        return Complex(self.re, self.im)

def deep_copy(p):
    #Copy of deepth 1 from a list
    copy = []
    for coef in p:
        copy.append(coef.copy())
    return copy


def mult_polynomial_complex(p, z):
    for coef in p:
        coef.mult_ip(z)

def add_polynomial(p1, p2):
    s = []
    m = min(len(p1), len(p2))
    d = max(len(p1), len(p2))
    if len(p1) < len(p2):
        p = p2
    else:
        p = p1

    for i in range(m):
        s.append(p1[i].add(p2[i]))
    for i in range(m, d):
        s.append(p[i].copy())
    return s


def developp(roots):
    p = [Complex(1,0)]
    for root in roots:
        temp = deep_copy(p)
        temp.insert(0, Complex(0,0))
        root.oppose()
        mult_polynomial_complex(p, root)
        root.oppose()

        p = add_polynomial(temp, p)
    return p

def print_polynomial(p):
    for i in range(len(p)):
        print("{} +i{} * x^{}".format(p[i].re, p[i].im, i))

def dist(z1, z2):
    w = Complex(z1.re - z2.re, z1.im - z2.im)
    return w.mod()

#print_polynomial(developp([Complex(-0.5, 1), Complex(10,0), Complex(5,0)]))

def is_free(z, l, eps):
    free = True 
    for w in l:
        if dist(w, z) <= eps:
            free = False
    return free

def random_complexe(maxmod = 1):
    z = Complex(random() * maxmod, random() * (2 * pi))
    while (abs(z.re) < 0.01) or (abs(z.im) < 0.01):
        z = Complex(random()* maxmod, random() * (2 * pi))
    return z


def create_test_ud(number, degree, eps):
    file_name = "D:/wamp/www/TIPE/Root_finder/test_set/test_set_ud/test_ud_{}_{}.txt".format(degree, eps)
    #eps refers here as the minimum distance beetween 2 roots
    
    file = open(file_name, "w")
    file.write(str(number) + "\n")
    for _ in range(number):
        roots = []
        i = 0
        while i < degree:
            z = random_complexe()
            if is_free(z, roots, eps) :
                i += 1
                roots.append(z)
        p = developp(roots)

        file.write(str(i) + "\n")
        for k in range(i):
            file.write("{:.15f}\n".format(p[k].re))
            file.write("{:.15f}\n".format(p[k].im))
    print(file_name)


def create_test_rd(number, degree):
    file_name = "D:/wamp/www/TIPE/Root_finder/test_set/test_set_rd/test_rd_{}.txt".format(degree)

    file = open(file_name, "w")
    file.write(str(number) + "\n")
    for _ in range(number):
        file.write(str(degree) + "\n")
        for _ in range(degree):
            z = random_complexe(maxmod=4)
            file.write("{:.15f}\n".format(z.re))
            file.write("{:.15f}\n".format(z.im))
    print(file_name)

def generate_test_ud():
    for i in range(10,30):
        create_test_ud(50, i, 1e-3)
#Degré très limité car au dela de 30, les coefficients sont petit à petit arondis à 0

def generate_test_rd():
    for i in range(20,151):
        create_test_rd(10, i)

#generate_test_ud()

def create_params(degmin, degmax, kmin, kmax, precision, EPS, name, mode, relaxation = 1):
    #EPS refers here to the asked precision, and precision to the minimal distance beetween 2 roots in case of mode = 'ud'

    filename = "D:/wamp/www/TIPE/Root_finder/" + name + ".txt"
    with open(filename, 'w', newline='', encoding='utf-8') as file:
        nb_test = (kmax - kmin + 1) * (degmax - degmin + 1)
        output_file = 'data/test_result_{0}/w_{1}/test_{0}_{2}.csv\n'.format(mode, relaxation, EPS)
        file.write("{}\n".format(nb_test))
        file.write(output_file)
        file.write("{}\n".format(EPS))
        file.write("{:.3f}\n".format(relaxation))
        for k in range(kmin, kmax + 1,20):
            for d in range(degmin, degmax+1):
                if mode == "ud":
                    test_file = 'test_set/test_set_ud/20/test_ud_{}_{}.txt\n'.format(d, precision)
                else:
                    test_file = 'test_set/test_set_rd/test_rd_{}.txt\n'.format(d)

                file.write(test_file)
                file.write("{}\n".format(k))
                

create_params(20,49,40,40,1e-03,1e-03, "input/test_rb_20_49_1e-03", "rb", 1.50)