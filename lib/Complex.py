from math import *

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
        
    def add_ip(self, z):
        if (type(z) != Complex):
            ValueError('Add parameter must be a Complex')
        else:
            self.re = self.re + z.re 
            self.im = self.im + z.im 
    
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

def horner(pol : list[Complex], z : Complex) -> Complex:
    if pol == []:
        return Complex(0,0)
    d = len(pol)-1
    eval = pol[d].copy()
    i = d-1
    while i >= 0:
        eval.mult_ip(z)
        eval.add_ip(pol[i])
        i -= 1
    return eval

#print_polynomial(developp([Complex(-0.5, 1), Complex(10,0), Complex(5,0)]))
