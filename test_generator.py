from random import*
from sys import*
from math import*

from Complex import*

'''file format V1:
    number of polynomials
    degree of P0
    Re(P0_1)
    IM(PO_1)
    Re(P0_2)
    .
    .
    .
    degree of P1
    Re(P1_1)
    ...
    The degree-th coef ins't shown, as it is equal to 1'''

'''file format V2:
    number of polynomials
    degree of P0
    Re(Root_1)+Im(Root_1)
    Re(Root_2)+Im(Root_2)
    .
    .
    .
    Re(P0_1)
    IM(PO_1)
    Re(P0_2)
    .
    .
    .
    degree of P1
    Re(P1_1)
    ...

    The degree-th coef ins't shown, as it is equal to 1'''

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


def create_test_ud_v1(number, degree, eps):
    file_name = "D:/wamp/www/TIPE/Root_finder/test_set/test_set_ud/test_v1_ud_{}_{}.txt".format(degree, eps)
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

def create_test_ud_v2(number, degree, eps):
    #eps refers here as the minimum distance beetween 2 roots
    file_name = "D:/wamp/www/TIPE/Root_finder/test_set/test_set_ud_v2/test_ud_v2_{}.txt".format(degree)
    
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
        for root in roots:
            file.write("{:.15f}+i{:.15f}\n".format(root.re, root.im))
        for k in range(i):
            file.write("{:.15f}\n".format(p[k].re))
            file.write("{:.15f}\n".format(p[k].im))
    print(file_name)


def create_test_rd_v1(number, degree):
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

def generate_test_ud_v1():
    for deg in range(10,30):
        create_test_ud_v1(50, deg, 1e-3)
#Degré très limité car au dela de 30, les coefficients sont petit à petit arondis à 0

def generate_test_ud_v2():
    for deg in range(10,11):
        create_test_ud_v2(15, deg, 1e-2)
#Degré très limité car au dela de 30, les coefficients sont petit à petit arondis à 0

def generate_test_rd():
    for deg in range(20,151):
        create_test_rd_v1(10, deg)


generate_test_ud_v2()


def create_params(version, degmin, degmax, kmin, kmax, precision, EPS, name, mode, relaxation = 1):
    #EPS refers here to the asked precision, and precision to the minimal distance beetween 2 roots in case of mode = 'ud'

    filename = "D:/wamp/www/TIPE/Root_finder/" + name + ".txt"
    with open(filename, 'w', newline='', encoding='utf-8') as file:
        nb_test = (kmax - kmin + 1) * (degmax - degmin + 1)
        if version == 1:
            output_file = 'data/test_result_{0}/w_{1}/test_{0}_{2}.csv\n'.format(mode, relaxation, EPS)
        if version == 2:
            output_file = 'data/test_result_{0}/test_{0}_v2_{1}.txt\n'.format(mode, EPS)
        file.write("{}\n".format(nb_test))
        file.write(output_file)
        file.write("{}\n".format(EPS))
        file.write("{:.3f}\n".format(relaxation))
        for k in range(kmin, kmax + 1,20):
            for d in range(degmin, degmax+1):
                if version == 1:
                    if mode == "ud":
                        test_file = 'test_set/test_set_ud/20/test_ud_{}_{}.txt\n'.format(d, precision)
                    else:
                        test_file = 'test_set/test_set_rd/test_rd_{}.txt\n'.format(d)
                if version == 2 and mode == 'ud':
                    test_file = 'test_set/test_set_ud_v2/test_ud_v2_{}.txt\n'.format(d)

                file.write(test_file)
                file.write("{}\n".format(k))
                

#create_params(2,10,10,40,40,1e-03,1e-03, "input/test_ud_v2_1e-03", "ud", 1.)
