# Euler axisymétrique — lot C

Solveur indépendant de SIMPLE, sans swirl, viscosité ni espèces. `Mesh` reste
un maillage méridien (x,r). `AxisymmetricEulerSolver` possède les métriques
`AxisymmetricGeometry`, une copie initiale des quatre tableaux SoA et les
buffers persistants nécessaires aux pas suivants. Le Mesh doit lui survivre.
`state()` expose uniquement l'état courant en lecture ; `step(t_final)` accepte
un pas et `advance_to(t_final)` atteint le temps demandé. Les vues d'état ne
doivent pas être conservées à travers un pas accepté.

L'assemblage est `V dU/dt = -sum(A Fhat) + G [0,0,p,0]`, avec les flux
Rusanov internes, le flux de paroi décrit ci-dessous, les
orientations owner/neighbor opposées et les métriques de révolution du lot A.
La source radiale compense le flux de pression au repos ; masse et énergie
n'ont pas de source. Les intégrales et les défauts de bilan de chaque pas
accepté sont disponibles dans `diagnostics()`, en distinguant flux sortants et
source radiale.

Chaque groupe doit recevoir une BC explicite : `Axis` uniquement sur r=0
(surface nulle, aucun flux), `SlipWall` par résolution de Riemann normale
(masse/énergie exactement imperméables), ou `Transmissive` par extrapolation intérieure.
Cette dernière n'est **pas** une condition caractéristique générale pour une
tuyère ; ici, les ondes de Sod n'atteignent pas les extrémités avant l'arrêt.

`compute_euler_slip_wall_flux()` résout le problème symétrique intérieur/miroir
avec une vitesse normale de paroi nulle, sans construire un état fantôme.
La normale est sortante : `un>0` comprime le gaz contre la paroi. Le flux est
exactement `[0, p_wall*nx, p_wall*nr, 0]`, indépendamment de la vitesse tangentielle.
La branche compression utilise Rankine–Hugoniot :
`c_shock = k*un + hypot(a,k*un)`, `k=(gamma+1)/4`,
`p_wall=p+rho*un*c_shock`. La branche expansion utilise les invariants
isentropiques : `p_wall=p*(1+(gamma-1)*un/(2*a))^(2*gamma/(gamma-1))`.
Ces relations sont celles des [courbes d'ondes Euler et du vide de
Riemann](https://www.clawpack.org/riemann_book/html/Euler.html).
Si `-un >= 2*a/(gamma-1)`, la solution exacte présente un vide à la paroi :
le flux et la pression y sont nuls, sans introduire de cellule vide ni de
plancher de pression. Une pression positive non représentable est au contraire
une erreur explicite, jamais transformée en vide.

La borne CFL de paroi inclut `abs(un)+a` (tête de raréfaction, y compris vide).
En compression, elle inclut aussi `c_shock` et `sqrt(gamma*p_wall/rho_inside)` :
la première borne la vitesse du choc dans le repère de la paroi ; la seconde
majore le son dans l'état comprimé car `rho_star >= rho_inside`.
La borne est arrondie vers le haut lors de sa conversion en double. Aucun
paramètre empirique ni allocation n'est ajouté aux calculs de face.

Euler explicite utilise un pas global `CFL min(V / sum(A lambda))`, CFL=0.2
par défaut. Cette estimation de faces ne garantit pas à elle seule
l'admissibilité avec la source géométrique : les candidats invalides entraînent
une division globale du pas par deux, avec un nombre de tentatives borné.
Aucun clipping n'est appliqué. Les erreurs de flux/assemblage ne sont pas
réessayées. Un échec conserve le dernier état accepté et signale cellule,
état candidat et pas ; une absence de progression temporelle est rejetée.

## Exemple Sod reproductible

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_VERIFICATION=ON
cmake --build build-release --target cfd_axisymmetric_euler_sod --parallel
./build-release/cfd_axisymmetric_euler_sod /tmp/euler-sod
```

Le driver résout le tube de longueur/rayon 1, gamma=1.4, R=1, avec les états
Sod `(rho,ux,p)=(1,0,1)` et `(0.125,0,0.1)`, séparés en x=0.5 ; ur=0.
Il compare N=100 et 200, quatre couches radiales QUAD, à t=0.2. La référence
indépendante utilise les invariants de raréfaction et Rankine–Hugoniot,
avec intersection des courbes d'ondes par dichotomie ([construction de
Riemann, Clawpack](https://www.clawpack.org/riemann_book/html/Euler.html)).
Les erreurs RMS/Linf sont évaluées aux centres des cellules et pondérées
par l'aire méridienne ; aucun ordre lisse n'est exigé aux discontinuités.
Les bilans intégrés utilisent, eux, les volumes de révolution.

Ouvrir `N100/sod.pvd` ou `N200/sod.pvd` dans ParaView. Chaque série contient
cinq VTU à t=0, 0.05, 0.1, 0.15, 0.2 : `rho`, `p`, `T`, `ux`, `ur`,
`velocity_magnitude`, `Mach`, `rhoE`. La représentation reste une coupe
méridienne, pas un maillage volumique 3D. Les champs dérivés ne sont calculés
qu'aux sorties, dans des buffers réutilisés.

Limites : ordre 1 diffusif, coût du pas explicite, précision énergétique
limitée aux très grands Mach (contrat du gaz parfait), pas de BC
caractéristiques ni de fonctionnement complet RPA. Les expansions radiales
et études axisymétriques avancées restent hors de ce lot.
