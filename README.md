# CFD Solver - starter

Squelette minimal d'un futur solveur CFD 2D en C++.

## Structure actuelle

- `Node` : point 2D.
- `Face` : arête reliant deux noeuds.
- `Cell` : type de cellule 2D. On commence avec les triangles, mais l'architecture prévoit aussi les quadrilatères.
- `Mesh` : conteneur principal du maillage.

Aucune génération de connectivité n'est encore implémentée.

## Compilation

Depuis la racine du projet :

```bash
cmake -S . -B build
cmake --build build
./build/cfd_solver
```

## Première étape de développement

Construire manuellement un maillage très simple composé de deux triangles :

```text
3 -------- 2
|        / |
|      /   |
|    /     |
|  /       |
0 -------- 1
```

Puis vérifier :
- les coordonnées des 4 noeuds ;
- la connectivité des 2 cellules ;
- ensuite seulement, construire automatiquement les faces.
