import RationalComplex as RC
import RCPolynomial as pol
import fractions as f
import math

def radius_complex_hull(p : pol.RCPolynomial) -> float:
    #Gershgorin (p must be monic)
    upper_bound = p.coefficients[0].mod()
    for i in range(1, p.degree):
        upper_bound = max(upper_bound, 1 + p.coefficients[i].mod())
    return max(1, upper_bound)

def simplify(p : pol.RCPolynomial) -> pol.RCPolynomial:
    normalizes_p = p.copy()
    normalizes_p.normalize()
    dp = normalizes_p.derivate()
    multiple_degree_zeros = pol.pgcd(normalizes_p, dp)
    simple_roots = pol.divide(normalizes_p, multiple_degree_zeros)
    return simple_roots

def Sd(degree : int, R : float, Njump : int) -> tuple[list[RC.RComplex], int]:
    nb_circles = math.ceil(0.26632 * math.log(degree))
    N = math.ceil(8.32547 * degree * math.log(degree))
    K = (1 + math.sqrt(2)) * R
    a = (degree - 1) / degree

    points_module = []
    point_argument = []

    for i in range(1, nb_circles + 1):
        points_module.append(K * math.pow(a, (2*i-1)/(4*nb_circles)))
    for k in range(Njump):
        for j in range(k,N,Njump):
            point_argument.append(2 * math.pi * j / N)

    Sd_set = [0 for _ in range(N*nb_circles)]
    for x in range(0,N):
        for y in range(0, nb_circles):
            re = points_module[y] * math.cos(point_argument[x])
            im = points_module[y] * math.sin(point_argument[x])
            Sd_set[y + x*nb_circles] = RC.RComplex(re, im)
        
    return Sd_set, N*nb_circles

def complex_zeroes(len : int) -> list[RC.RComplex]:
    tab = [RC.RComplex(0,0) for _ in range(len)]
    return tab 

def print_complex_array(arr : list[RC.RComplex]):
    print("[", end='')
    for i in range(len(arr)-1):
        arr[i].print()
        print(", ", end='')
    arr[len(arr)-1].print()
    print("]")

def add_to_set(roots : list[RC.RComplex], z : RC.RComplex, EPS : float) -> bool:
    new_root = True 
    for root in roots:
        if RC.distance(root, z) <= EPS:
            new_root = False
    if new_root:
        roots.append(z)
        return True
    else:
        return False

def step(z0 : RC.RComplex, prev : RC.RComplex, p : pol.RCPolynomial, dp : pol.RCPolynomial, w : float):
    prev.replace(z0)

    pz = p.evaluate(z0)
    dpz = dp.evaluate(z0)
    try:
        pz.divide_ip(dpz)
    except ZeroDivisionError:
        raise RC.Infinity
    pz.mult_scalar_ip(w)
    pz.opposite_ip()
    pz.add_ip(z0)

    z0.replace(RC.RComplex(f.Fraction(float(pz.re), f.Fraction(float(pz.im)))))

def evaluate_mod(p : pol.RCPolynomial, z : RC.RComplex) -> float:
    pz = p.evaluate(z)
    return pz.mod()

def solve_deg2(p : pol.RCPolynomial) -> list[RC.RComplex]:
    '''
    p must pe monic
    '''
    if p.degree == 2:
        b2 = RC.mult(p.coefficients[1], p.coefficients[1])
        delta2 = p.coefficients[0].copy()
        delta2.mult_ip(p.coefficients[2])
        delta2.mult_scalar_ip(4.)
        delta2 = RC.sub(b2, delta2)
        delta2.opposite_ip()
        delta = RC.sqrt(delta2)

        root0 = p.coefficients[1].opposite()
        root0.sub_ip(delta)
        root0.mult_scalar_ip(0.5)

        root1 = p.coefficients[1].opposite()
        root1.add_ip(delta)
        root1.mult_scalar_ip(0.5)

        return [root0, root1]
    elif p.degree == 1:
        root = p.coefficients[0].copy()
        root.divide_ip(p.coefficients[1])

        return [root]
    else:
        raise ValueError


# def main():
#     roots = []
#     add_to_set(roots, RC.RComplex(1,0), 1.)
#     print_complex_array(roots)

#     S2, size = Sd(2,1.,1)
#     print_complex_array(S2)

#     roots = [RC.RComplex(1,0), RC.RComplex(0.5,0)]
#     p = pol.developp(roots, 2)
#     p.print()
#     print("")

#     print_complex_array(solve_deg2(p))

# main()
