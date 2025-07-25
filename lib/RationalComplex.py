import math
import fractions as f

class Infinity(Exception):
    pass

class RComplex:
    def __init__(self,x,y):
        self.re = f.Fraction(x).limit_denominator()
        self.im = f.Fraction(y).limit_denominator()

    def print(self):
        print(self.re.__format__(".16f"), end='')
        print(" +i", end='')
        print(self.im.__format__(".16f"), end='')

    def copy(self) -> "RComplex":
        return RComplex(self.re, self.im)
    
    def reset(self):
        self.re = f.Fraction()
        self.im = f.Fraction()

    def replace(self, z : "RComplex"):
        self.re = f.Fraction(z.re)
        self.im = f.Fraction(z.im)

    def mod(self) -> float:
        m = self.re * self.re + self.im * self.im
        return math.sqrt(m)
    
    def arg(self) -> float:
        if self.re == 0:
            if self.im > 0:
                return math.pi / 2
            else:
                return - math.pi / 2
        else:
            return math.atan(float(self.im / self.re))
    
    def mod2(self) -> f.Fraction:
        m = self.re * self.re + self.im * self.im
        return m.limit_denominator()
    
    def conjugate(self) -> "RComplex":
        return RComplex(self.re, - self.im)
    
    def conjugate_ip(self):
        self.im = - self.im

    def opposite(self) -> "RComplex":
        return RComplex(- self.re, - self.im)
    
    def opposite_ip(self):
        self.re = - self.re
        self.im = - self.im

    def add_ip(self, z : "RComplex"):
        self.re = self.re + z.re
        self.im = self.im + z.im

    def sub_ip(self, z : "RComplex"):
        self.re = self.re - z.re
        self.im = self.im - z.im

    def mult_ip(self, z : "RComplex"):
        RE = self.re * z.re - self.im * z.im
        IM = self.re * z.im + self.im * z.re
        self.re = RE
        self.im = IM

    def mult_scalar_ip(self, x):
        self.re = self.re * f.Fraction(x).limit_denominator()
        self.im = self.im * f.Fraction(x).limit_denominator()

    def invert(self) -> "RComplex":
        m = self.re * self.re + self.im * self.im
        if m == 0:
            raise ZeroDivisionError
        return RComplex(self.re / m, - self.im / m)

    def invert_ip(self):
        m = self.re * self.re + self.im * self.im
        if m == 0:
            raise ZeroDivisionError
        self.re = self.re / m
        self.im = - self.im / m

    def divide_ip(self, z : "RComplex"):
        m = self.re * self.re + self.im * self.im
        if m == 0:
            raise ZeroDivisionError
        temp = z.conjugate()
        RE = self.re * temp.re - self.im * temp.im
        IM = self.re * temp.im + self.im * temp.re
        self.re = RE / m
        self.im = IM / m

    def is_null(self):
        return (self.re == 0 and self.im == 0)
    
def conjugate(z : RComplex):
    return RComplex(z.re, -z.im)

def opposite(z : RComplex):
    return RComplex(- z.re, - z.im)

def invert(z : RComplex):
    m = z.re * z.re + z.im * z.im
    if m == 0:
        raise ZeroDivisionError
    return RComplex(z.re / m, - z.im / m)

def add(z1 : RComplex, z2 : RComplex) -> RComplex:
    return RComplex(z1.re + z2.re, z1.im + z2.im)

def sub(z1 : RComplex, z2 : RComplex) -> RComplex:
    return RComplex(z1.re - z2.re, z1.im - z2.im)

def mult(z1 : RComplex, z2 : RComplex) -> RComplex:
    RE = z1.re * z2.re - z1.im * z2.im
    IM = z1.re * z2.im + z1.im * z2.re
    return RComplex(RE, IM)

def divide(z1 : RComplex, z2 : RComplex) -> RComplex:
    m = z2.re * z2.re + z2.im * z2.im
    if m == 0:
        raise ZeroDivisionError
    temp = z2.conjugate()
    RE = z1.re * temp.re - z1.im * temp.im
    IM = z1.re * temp.im + z1.im * temp.re
    return RComplex(RE / m, IM / m)

def sqrt(z : RComplex) -> RComplex:
    r = math.sqrt(z.mod())
    arg = z.arg() / 2
    return RComplex(r * math.cos(arg), r * math.sin(arg))

def distance(z1 : RComplex,z2 : RComplex) -> float:
    temp = sub(z2, z1)
    return temp.mod()

# def main():
#     z1 = RComplex(0.25,-0.5)
#     z1.print()
#     print("\n", end="")
#     z2 = RComplex(1,2)
#     z1.conjugate_ip()
#     z1.add_ip(z2)
#     z1.print()
#     print("\n", end="")
#     z3 = divide(z1,z2)
#     z3.print()
#     print("\n", end="")
# main()
