import numpy as np
import matplotlib.pyplot as plt
from math import *
from collections import defaultdict
import csv
from Complex import *

root_path = ""

def unique_value(attribute, data):
    values = []
    for test in data:
        values.append(test[attribute])
    return np.unique(sorted(values))

def attribute_mean(attribute, data):
    s = 0
    nb_value = 0
    for test in data:
        nb_value +=1
        s += float(test[attribute])
    return s / nb_value

def attribute_mean_by_group(attribute, data, group, group_name='degree'):
    s = dict(zip(group, [0 for k in range(len(group))]))
    nb = dict(zip(group, [0 for k in range(len(group))]))
    for test in data:
        s[int(test[group_name])] += float(test[attribute])
        nb[int(test[group_name])] += 1
    res = []
    keys = sorted(s.keys()) 
    for key in keys:
        res.append(s[key] / nb[key])
    return np.array(res)

def get_attribute(attribute, data):
    values = []
    for test in data:
        values.append(float(test[attribute]))
    return values

def fetch(files, version, precision):
    nb_failed = 0
    nb_tot = 0
    succes = []
    if version == 1:
        for file in files:
            with open(file, "r", newline='', encoding='utf-8') as f:
                data = csv.DictReader(f, ["k", "max_error", "degree", "NM_iterations", "exec_time"], delimiter=';')
                for test in data:
                    if test['max_error'] == 'inf' or(float(test['max_error']) > precision or float(test['max_error']) < 0):
                        nb_failed +=1
                    else:
                        succes.append(test)
                    nb_tot += 1
        succes_rate = (nb_tot - nb_failed) / nb_tot
        return succes, succes_rate
    if version == 2:
        for file in files:
            f = open(file, "r", encoding='utf-8')
            line = f.readline()
            while True:
                if line == "New test file\n":
                    test_filename = root_path + f.readline().strip()
                    test_file = open(test_filename, "r", encoding='utf-8')
                    nb_test = int(test_file.readline())
                    for i in range(nb_test):
                        #get data
                        line = f.readline().strip() 
                        nb_tot += 1

                        #if test failed
                        if line == "-1":
                            nb_failed += 1
                            d = int(test_file.readline().strip())
                            for _ in range(3*d):
                                test_file.readline()
                            continue

                        #if test didn't failed
                        perf = line.split(";")
                        test = {"k" : int(perf[0]),
                                "degree" : int(perf[1]),
                                "NM_iterations" : int(perf[2]),
                                "exec_time" : float(perf[3]),
                                "max_error" : -1}
                        roots_founded = []
                        for _ in range(test["degree"]):
                            root = f.readline().strip().split("+i")
                            roots_founded.append(Complex(float(root[0]), float(root[1])))

                        #get polynomial
                        coefficients = []
                        d = int(test_file.readline().strip())
                        if d != test['degree']:
                            raise ValueError
                        for _ in range(d):
                            test_file.readline()
                        for _ in range(d):
                            re = float(test_file.readline().strip())
                            im = float(test_file.readline().strip())
                            coefficients.append(Complex(re,im))
                        coefficients.append(Complex(1,0))

                        #compute max_error
                        for approx_root in roots_founded:
                            z = horner(coefficients, approx_root)
                            if z.mod() >= test["max_error"]:
                                test["max_error"] = z.mod()
                            
                        #add test
                        succes.append(test)
                else:
                    break
            f.close()

        succes_rate = (nb_tot - nb_failed) / nb_tot
        return succes, succes_rate
                

def succes_rate(files, version, precision):
    _, succes_rate = fetch(files, version, precision)
    print("Succes rate : {:.3f}".format(succes_rate))

def chart(files, version, precision):
    data, _ = fetch(files, version, precision)
    
    degree = unique_value('degree', data)
    degree = np.sort(np.array(list(map(int, degree))))
    mean_exec_time = attribute_mean_by_group('exec_time', data, degree)
    mean_nm_it = attribute_mean_by_group('NM_iterations', data, degree)
    exec_time = get_attribute('exec_time', data)
    nb_iteration = get_attribute('NM_iterations', data)

    group = defaultdict(lambda: defaultdict(list))
    for row in data:
        try:
            k = int(row["k"])
            nm_iterations_k = float(row["NM_iterations"])
            exec_time_k = float(row['exec_time'])
            group[k][0].append(nm_iterations_k)
            group[k][1].append(exec_time_k)
        except (ValueError, KeyError):
            continue

    fig, axes = plt.subplots(2,2)
    axes[0, 0].plot(degree, mean_exec_time, label="Temps d'exécution moyen")
    axes[0, 0].set_xlabel("Degré du polynôme")
    axes[0, 0].set_ylabel("Temps d'exécution moyen (s)")
    axes[0, 0].set_title("Temps d'exécution moyen en fonction du degré du polynôme")
    axes[0, 0].grid()
    axes[0, 0].legend()

    axes[1, 0].plot(degree, mean_nm_it, label="Itérations moyennes")
    axes[1, 0].set_xlabel("Degré du polynôme")
    axes[1, 0].set_ylabel("Nombre d'itérations moyen")
    axes[1, 0].set_title("Nombre d'itérations moyen de la méthode en fonction du degré du polynôme")
    axes[1, 0].grid()
    axes[1, 0].legend()

    axes[0, 1].scatter(exec_time, nb_iteration, label="Itérations")
    axes[0, 1].set_xlabel("Temps d'exécution (s)")
    axes[0, 1].set_ylabel("Nombre d'itérations")
    axes[0, 1].set_title("Nombre d'itérations en fonction du temps d'exécution")
    axes[0, 1].grid()
    axes[0, 1].legend()

    for k in group.keys():
        axes[1, 1].scatter(np.mean(group[k][0]), np.mean(group[k][1]), marker='o', label=f'k={k}')
    axes[1, 1].set_xlabel("Average NM_iterations")
    axes[1, 1].set_ylabel("Average exec_time")
    axes[1, 1].set_title("Average exec_time vs Average NM_iterations for each k")
    axes[1, 1].grid()

    plt.tight_layout()
    plt.show()

def exec_time_vs_nmi_vs_k(files, version, precision):
    succes, _ = fetch(files, version, precision)
    group = defaultdict(lambda: defaultdict(list))
    
    for row in succes:
        try:
            k = int(row["k"])
            nm_iterations = float(row["NM_iterations"])
            exec_time = float(row['exec_time'])
            group[k][0].append(nm_iterations)
            group[k][1].append(exec_time)
        except (ValueError, KeyError):
            continue
    
    for k in group.keys():        
        plt.scatter(np.mean(group[k][0]), np.mean(group[k][1]), marker='o', label=f'k={k}')
    
    plt.xlabel("Average NM_iterations")
    plt.ylabel("Average exec_time")
    plt.title("Average exec_time vs Average NM_iterations for each k")
    plt.legend()
    plt.grid()
    plt.show()

def degree_vs_nmi_vs_k(files, version, precision):
    succes, _ = fetch(files, version, precision)
    group = defaultdict(lambda: defaultdict(list))
    
    for row in succes:
        try:
            k = int(row["k"])
            degree = int(row["degree"])
            nm_iterations = float(row["NM_iterations"])
            group[k][degree].append(nm_iterations)
        except (ValueError, KeyError):
            continue
    
    for k, degree_data in group.items():
        degrees = sorted(degree_data.keys())
        means = [np.mean(degree_data[d]) for d in degrees]
        
        plt.plot(degrees, means, marker='o', label=f'k={k}')
    
    plt.xlabel("Degree")
    plt.ylabel("Average NM_iterations")
    plt.title("Average NM_iterations vs Degree for each k")
    plt.legend()
    plt.grid()
    plt.show()

def degree_vs_exec_time_vs_k(files, version, precision):
    succes, _ = fetch(files, version, precision)
    group = defaultdict(lambda: defaultdict(list))
    
    for row in succes:
        try:
            k = int(row["k"])
            degree = int(row["degree"])
            exec_time = float(row["exec_time"])
            group[k][degree].append(exec_time)
        except (ValueError, KeyError):
            continue
    
    for k, degree_data in group.items():
        degrees = sorted(degree_data.keys())
        means = [np.mean(degree_data[d]) for d in degrees]
        
        plt.plot(degrees, means, marker='o', label=f'k={k}')
    
    plt.xlabel("Degree")
    plt.ylabel("Average exec_time")
    plt.title("Average exec_time vs Degree for each k")
    plt.legend()
    plt.grid()
    plt.show()

def degree_vs_exec_time_vs_file(files, version, precision):
    i = 0
    for file,label in files:
        data, _ = fetch([file], version, precision)
    
        degrees = unique_value('degree', data)
        degrees = np.sort(np.array(list(map(int, degrees))))
        mean_exec_time = attribute_mean_by_group('exec_time', data, degrees)
            
        plt.plot(degrees, mean_exec_time, label=f'{label}')
        i+=1
    
    plt.xlabel("Degree")
    plt.ylabel("Temps d'exécution (s)")
    plt.title("Temps d'exécution en fonction du degré du polynôme, par valeur de w")
    plt.legend()
    plt.grid()
    plt.show()


succes_rate([root_path + "/data/test_result_ud/test_ud_v2_0.001.txt"], 2, 1e-03)
#chart([root_path + "/data/Njump/test_ud_1.csv"], 1, 0.001)
#chart([
#     root_path + "/data/test_ud_eval_0.001_cst.csv",
#     root_path + "/data/test_ud_eval_1e-06_cst.csv",
#     root_path + "/data/test_ud_eval_1e-09_cst.csv"], 1, 0.001)
#degree_vs_exec_time_vs_file([(root_path + "/data/w/test_ud_w_{}.csv".format(k),k) for k in [0.25,0.32,0.5,0.65,0.85,1.15,1.3,1.5]], 0.001)
#degree_vs_exec_time_vs_k([root_path + "/data/test_ud_20_1e-09.csv"], 1e-09)

def tchebychev(n):
    dmin = 1
    for k in range(2, ceil(n/2)):
        d = abs(cos((2*k-1)*pi/(2*n)) - cos((2*k-3)*pi/(2*n)))
        if d < dmin:
            dmin = d
    print(dmin)

# tchebychev(50)
# tchebychev(100)
# tchebychev(150)
def find_nearest(x, list):
    nearest = list[0]
    for elem in list:
        if abs(x - elem) < abs(x - nearest):
            nearest = elem
    return nearest
    
def result_test_tchebychev(file):
    f =  open(file, "r", newline='', encoding='utf-8')
    for _ in range(14):
        info = (f.readline()).split(";")
        n = int(info[0])
        nb_root = int(info[1])
        exec_time = float(info[2])
        roots = []
        if n == nb_root:
            for _ in range(n):
                roots.append(float(f.readline()))
            values = np.array([cos((2*k-1)*pi/(2*n)) for k in range(1, n+1)])
            associate = np.array([find_nearest(val, roots) for val in values])
            d = np.max(np.abs(values - associate))
            print("Erreur maximale pour T{} : {}.\nTemps d'execution : {} s.\n".format(n, d, exec_time))

#result_test_tchebychev(root_path + "/data/tchebychev.csv")

def result_rb(files, precision):
    degree = np.array([d for d in range(20, 50)])

    nb_failed = {k:0 for k in degree}
    nb_tot = {k:0 for k in degree}
    succes = {d:np.array([]) for d in degree}
    for file in files:
        with open(file, "r", newline='', encoding='utf-8') as f:
            data = csv.DictReader(f, ["k", "max_error", "degree", "NM_iterations", "exec_time"], delimiter=';')
            for test in data:
                if test['max_error'] == 'inf' or(float(test['max_error']) > precision or float(test['max_error']) < 0):
                    nb_failed[int(test['degree'])] +=1
                else:
                    succes[int(test['degree'])] = np.append(succes[int(test['degree'])], float(test['exec_time']))
                nb_tot[int(test['degree'])] +=1
    
    mean_exec_time = []
    for k in degree:
        if succes[k].size > 0:
             mean_exec_time.append(np.mean(succes[k]))
        else:
            mean_exec_time.append(0)
    succes_rate = [1 - nb_failed[k]/nb_tot[k] for k in degree]

    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.bar(degree, succes_rate)
    ax1.set_xlabel("Degré du polynôme")
    ax1.set_ylabel("Taux de réussite")
    ax1.set_title("Taux de réussite par degré")
    ax1.grid()

    ax2.bar(degree, mean_exec_time)
    ax2.set_xlabel("Degré du polynôme")
    ax2.set_ylabel("Temps d'exécution moyen (s)")
    ax2.set_title("Temps d'exécution moyen en fonction du degré du polynôme")
    ax2.grid()

    plt.tight_layout()
    plt.show()

#result_rb([root_path + "/data/test_rb_20_49_0.001.csv"], 0.001)
