# PolynomialRootFinder
Algorithme de recherche de racines d'un polynôme utilisant les résultats démontrés [ici](https://www.math.stonybrook.edu/~scott/Papers/Newton-HSS.pdf).

# Implémentation
Les programmes sont écrits en C et en Python.
On dispose ici :
- d'un générateur de fichiers tests
- de librairie auxiliaires permettant de pouvoir effectuer des calculs avec des complexes et des polynôme de $\mathbb{C}[X]$, de lire les fichiers de test
- un fichier __aux.c__ pour les fonctions auxiliares du programme : le calcul de $S_d$ par exemple
- Plusieurs implémentation de l'algorithme se différenciant par leur condition d'arrêt, la possibilité ou non de les utiliser sur des polynômes dont les racines ne sont pas contenues dans le disque unité
- Une implémentation spécifique pour les polynômes de Tchebychev
- Un fichier permettant d'analyser les résultats obtenus


L'implémentation est la suivante :
- On calcule $\frac{P}{pgdc(P,P')}$ pour avoir un polynôme unitaire scindé à racine simple
- On calcule $S_d$
- On applique la méthode de Newton successivement sur chaque éléments de $S_d$, en passant à l'élément suivant après au plus $K$ itérations
- Si on a parcouru tout l'ensemble $S_d$ sans avoir trouvé toutes les racines de $p$, on recommence en effectuant de nouveau au plus $K$ itérations sur chaque point de $S_d$, et ainsi de suite ...

# Améliorations à venir
Les tests montrent que les calculs avec des nombres à virgule flottante sont beaucoup trop imprécis dans se contexte (ou l'évaluation du polynôme doit être très faible, donc il faut que les termes se compensent).

Pour mieux tester l'algorithme, une implémentation Python calculant des approximations des complexes de l'ensemble $S_d$ et des coefficients des polynômes par des complexes à parties réelles et imaginaires rationnelles, et n'effectuant que des calculs exactes est à venir.
