import fractions as f
import RationalComplex as RC

def indice(n):
    return range(0,n+1)

class RCPolynomial:
    def __init__(self, d : int):
        self.degree = d
        self.coefficients = [RC.RComplex(0,0) for _ in indice(d)]
    
    def get(self, i : int) ->  RC.RComplex:
        return self.coefficients[i].copy()
    
    def set(self, i : int, z : RC.RComplex):
        self.coefficients[i] = z.copy()

    def print(self):
        if self.degree == -1:
            print("0")
        else :
            self.coefficients[0].print()
            for i in range(1, self.degree+1):
                print(" + ", end='')
                self.coefficients[i].print()
                print("X^{}".format(i), end='')

    def reset(self):
        self.degree = -1
        self.coefficients = []

    def copy(self):
        p = RCPolynomial(self.degree)
        for i in indice(self.degree):
            p.set(i, self.coefficients[i])
        return p
    
    def replace(self, p : "RCPolynomial"):
        self.degree = p.degree
        self.coefficients = [RC.RComplex(0,0) for _ in indice(p.degree)]
        for i in indice(p.degree):
            self.set(i, p.coefficients[i])

    def adjust_degree(self):
        while self.degree >=0 and self.coefficients[self.degree].is_null():
            self.degree = self.degree - 1
            self.coefficients.pop()
        
    def evaluate(self, z : RC.RComplex) -> RC.RComplex:
        if self.degree == -1:
            return RC.RComplex(0,0)
        res = self.coefficients[self.degree].copy()
        for i in range(self.degree - 1, -1, -1):
            res.mult_ip(z)
            res.add_ip(self.coefficients[i])
        return res

    def derivate(self):
        if self.degree < 1:
            return RCPolynomial(-1)
        dp = RCPolynomial(self.degree - 1)
        for i in indice(self.degree - 1):
            z = self.coefficients[i+1].copy()
            z.mult_scalar_ip(i+1)
            dp.set(i, z)
        return dp

    def mult_by_complex(self, z : RC.RComplex):
        if z.is_null():
            self.reset()
        else:
            for i in indice(self.degree):
                (self.coefficients[i]).mult_ip(z)

    def mult_by_scalar(self, x : f.Fraction):
        if x == 0:
            self.reset()
        else:
            for i in indice(self.degree):
                self.coefficients[i].mult_scalar_ip(x)

    def mult_X(self):
        if self.degree != -1:
            self.coefficients.insert(0, RC.RComplex(0,0))
            self.degree = self.degree + 1

    def normalize(self):
        if self.degree != -1:
            self.mult_by_complex(self.coefficients[self.degree].invert())

def add(p1 : RCPolynomial, p2 : RCPolynomial) -> RCPolynomial:
    if p1.degree == -1:
        return p1.copy()
    if p2.degree == -1:
        return p2.copy()
    
    if p1.degree <= p2.degree:
        small = p1
        big = p2 
    else:
        small = p2
        big = p1
    sum = RCPolynomial(big.degree)
    for k in indice(small.degree):
        sum.set(k, RC.add(small.coefficients[k], big.coefficients[k]))
    for k in range(small.degree+1, big.degree+1):
        sum.set(k, big.coefficients[k])
    return sum

def sub(p1 : RCPolynomial, p2 : RCPolynomial) -> RCPolynomial:
    if p1.degree == -1:
        diff = p1.copy()
        for k in indice(p1.degree):
            diff.set(k, diff.coefficients[k].opposite())
        return diff
    if p2.degree == -1:
        return p2.copy()
    
    if p1.degree <= p2.degree:
        small = p1
        big = p2 
    else:
        small = p2
        big = p1
    sum = RCPolynomial(big.degree)
    for k in indice(small.degree):
        sum.set(k, RC.sub(small.coefficients[k], big.coefficients[k]))
    for k in range(small.degree+1, big.degree+1):
        sum.set(k, big.coefficients[k].opposite())
    sum.adjust_degree()
    return sum

def mult(p1 : RCPolynomial, p2 : RCPolynomial) -> RCPolynomial:
    if p1.degree == -1 or p2.degree == -1:
        return RCPolynomial(-1)
    prod = RCPolynomial(p1.degree, p2.degree)
    for k in indice(p1.degree):
        for l in indice(p2.degree):
            prod.set(k+l, RC.add(prod.coefficients[k+l], RC.mult(p1.coefficients[k], p2.coefficients[l])))
    return prod

def modulus(a : RCPolynomial, b : RCPolynomial) -> RCPolynomial:
    if b.degree == -1:
        raise ZeroDivisionError
    
    p = a.copy()
    while p.degree >= b.degree:
        z = RC.divide(p.coefficients[p.degree], b.coefficients[b.degree])

        for i in indice(b.degree):
            t = RC.mult(b.coefficients[b.degree- i], z)
            p.set(p.degree - i, RC.sub(p.coefficients[p.degree- i], t))
        p.adjust_degree()
    return p

def pgcd(a : RCPolynomial, b : RCPolynomial) -> RCPolynomial:
    a_copy = a.copy()
    b_copy = b.copy()

    while b_copy.degree != -1:
        r = modulus(a_copy, b_copy)
        a_copy.replace(b_copy)
        b_copy.replace(r)
    a_copy.normalize()
    return a_copy

def divide(a : RCPolynomial, b : RCPolynomial) -> RCPolynomial:
    if b.degree == -1:
        raise ZeroDivisionError

    if b.degree == 0:
        q = a.copy()
        q.mult_by_complex(b.coefficients[0].invert())

    if a.degree < b.degree:
        return RCPolynomial(-1)
    
    r = a.copy()
    q = RCPolynomial(a.degree - b.degree)
    while r.degree >= b.degree:
        fac = RC.divide(r.coefficients[r.degree], b.coefficients[b.degree])
        q.set(r.degree - b.degree, fac)

        for i in indice(b.degree):
            t = RC.mult(b.coefficients[b.degree - i], fac)
            t = RC.sub(r.coefficients[r.degree - i], t)
            r.set(r.degree - i, t)
        r.adjust_degree()
    return q

def developp(roots : list[RC.RComplex], nb_roots : int) -> RCPolynomial:
    p = RCPolynomial(0)
    p.set(0, RC.RComplex(1,0))

    for k in range(0, nb_roots):
        temp = p.copy()
        p.mult_X()
        temp.mult_by_complex(roots[k])
        p = sub(p, temp)
    return p

# def main():
#     p = RCPolynomial(0)
#     p.print()
#     print("\n", end="")
#     p.set(0, RC.RComplex(1,0))
#     p.print()
#     print("\n", end="")
#     p.mult_X()
#     p.mult_by_complex(RC.RComplex(2,0))
#     p.print()
#     print("\n", end="")
#     roots = [RC.RComplex(1,0), RC.RComplex(0.5,0)]
#     p = developp(roots, 2)
#     p.print()
#     print("\n", end="")
#     dp = p.derivate()
#     dp.print()
#     print("\n", end="")
#     dp.normalize()
#     dp.print()
#     print("\n", end="")
#     t1 = pgcd(p,dp)
#     t1.print()
#     print("\n", end="")
#     t2 = divide(p,dp)
#     t2.print()
#     print("\n", end="")

# main()