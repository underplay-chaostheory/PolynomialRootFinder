import lib/RationalComplex as RC
import lib/RCPolynomial as pol
import lib/RootFinderLib as lib
import math

def root_finder(p : pol.RCPolynomial, EPS : float, ITERATION_LOOP : int, w : float = 1.0, N_jump : int = 1) -> list[RC.RComplex]:
    '''
    The roots of p must at least 3*EPS far of each others
    '''
    simple_roots = lib.simplify(p)
    dsr = simple_roots.derivate()
    if p.degree < 3:
        return lib.solve_deg2(simple_roots)
    
    nb_root = simple_roots.degree
    nb_root_found = 0
    roots = []

    upper_bound = lib.radius_complex_hull(simple_roots)
    Sd, size = lib.Sd(nb_root, upper_bound, N_jump)
    # lib.print_complex_array(Sd)
    set = Sd 
    aux = []
    set_size = size 
    aux_size = 0

    completed = 0
    nm_iteration = 0
    max_nm_iteration = set_size * 200
    tol = EPS / nb_root
    while nb_root_found < nb_root and completed < size and nm_iteration < max_nm_iteration:
        i = 0
        while i < set_size and nb_root_found < nb_root and completed < size:
            prev = RC.RComplex(0,0)
            z0 = set[i].copy()

            try:
                lib.step(z0,prev,simple_roots,dsr,w)
                k = 1
                nm_iteration += 1
                while k < ITERATION_LOOP and RC.distance(z0,prev) > tol:
                    lib.step(z0,prev,simple_roots,dsr,w)

                    k += 1
                    nm_iteration += 1
                
                if nb_root_found < nb_root and RC.distance(z0,prev) <= tol:
                    completed += 1
                    if lib.add_to_set(roots, prev, 2 * EPS):
                        nb_root_found += 1
                        prev.print()
                else:
                    aux.append(z0)
                    aux_size += 1
            except RC.Infinity:
                completed += 1
            
            i += 1
        
        set_size = aux_size
        aux_size = 0 
        set = aux 
        aux = []
    
    return roots

def main():
    p = pol.RCPolynomial(5)
    p.set(0, RC.RComplex(-1,0))
    p.set(5, RC.RComplex(1,0))

    roots = root_finder(p, 0.00001, 40, 1., 1)

    lib.print_complex_array(roots)

main()
