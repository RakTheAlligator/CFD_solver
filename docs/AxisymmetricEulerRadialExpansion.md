# Expansion radiale axisymétrique — lot D

## Solution exacte et équations

Posons `s=1+alpha0*t`, `a=alpha0/s`, `Q=p/(gamma-1)` et
`K=rho*a²*r²/2`. Avec `ux=0`, `ur=a*r`, `rho=rho0/s²` et
`p=p0/s^(2*gamma)`, on a `a'=-a²`, `rho'=-2*a*rho`,
`p'=-2*gamma*a*p` et `div(u)=2*a`.

- Continuité : `rho' + (1/r)*d(r*rho*ur)/dr = -2*a*rho+2*a*rho=0`.
- Momentum axial : tous les termes sont nuls car `ux=0` et `dp/dx=0`.
- Momentum radial : `d(ur)/dt+ur*d(ur)/dr=-a²*r+a²*r=0`, et `dp/dr=0`.
  En forme conservative, `d(rho*ur)/dt=-3*rho*a²*r` et
  `(1/r)*d[r*(rho*ur²+p)]/dr=3*rho*a²*r+p/r` ; le terme `p/r`
  du membre droit compense donc exactement la partie géométrique de pression.
- Énergie : `d(Q+K)/dt=-2*gamma*a*Q-4*a*K` et
  `(1/r)*d[r*(Q+K+p)*ur]/dr=2*gamma*a*Q+4*a*K`.
  Il n'existe aucune source physique de masse ou d'énergie totale.

Ces expressions sont régulières sur l'axe en forme intégrée : le flux à
`r=0` est nul. L'intégration du terme radial donne `2*pi*A_meridienne*p`,
conforme au coefficient `G` du lot A. Le mouvement radial et les variations
temporelles de densité/pression sont non triviaux, malgré l'accélération
matérielle nulle. Aucune équation source supplémentaire n'est introduite.

## Projection conservative, indépendante du solveur

Pour un rectangle de largeur axiale `dx`, `ri<=r<=ro`, le poids est
`dV=2*pi*r*dx*dr`. Posons `h=ro-ri` et `m=ri+h/2`. Alors

```
V = pi*dx*h*(ro+ri)
mu1 = <r>_V  = m + (h/m)*(h/12)
mu2 = <r²>_V = m² + h²/4
Var(r) = h²/12 * [1-(h/m)²/12]
```

Les différences de puissances presque égales sont ainsi évitées, y compris
pour les anneaux minces. À l'axe : `mu1=2*ro/3`, `mu2=ro²/2`.
La référence des tableaux SoA est exactement

```
Ubar = [rho, 0, rho*a*mu1, p/(gamma-1)+rho*a²*mu2/2].
```

La pression de l'EOS appliquée à ces moyennes est
`p_EOS(Ubar)=p+(gamma-1)*rho*a²*Var(r)/2`, et non la pression
ponctuelle uniforme `p`. Les deux comparaisons seront rapportées séparément.
Une quadrature de Gauss indépendante vérifiera les moments et la projection.

## Frontières et protocole

Le premier domaine étudié est `[0,1]x[0,1]`, avec Axis en bas et Transmissive sur
les trois autres côtés. L'invariance axiale est compatible avec l'extrapolation
aux extrémités. Une SlipWall extérieure serait incompatible avec `ur(R)>0`.

Transmissive utilise la dernière moyenne cellulaire, pas la trace exacte en R.
Au temps initial, son flux de masse extérieur relatif à l'exact vaut
`mu1_last/R` : son défaut relatif est donc `1-mu1_last/R=O(h/R)`.
L'erreur d'enthalpie et la viscosité Rusanov contribuent également à l'erreur.
Les normes principales incluent **toutes** les cellules, y compris l'axe et
la frontière extérieure. Aucun retrait de cellules ne sert à faire passer
la vérification. Un contrôle sur un domaine étendu, avec même maille et
mêmes pas dans la partie commune, mesure la sensibilité à la frontière.

Paramètres principaux : `rho0=p0=alpha0=1`, `gamma=1.4`, `t_final=0.2`,
`Nr=20,40,80`, `Nx=Nr/4`. Le raffinement porte sur le même domaine physique.
Les normes sont `L1=sum(V*abs(e))/sum(V)`,
`RMS=sqrt(sum(V*e²)/sum(V))`, `Linf=max(abs(e))`, pour rho, rho*ur et rhoE.
L'ordre observé est `log2(E_coarse/E_fine)` pour un doublement de Nr.
À CFL fixé, dt se raffine également avec h : il s'agit d'une convergence
couplée espace/temps, pas d'une isolation de l'erreur spatiale pure.
Sur le maillage fixé, les CFL 0.4, 0.2, 0.1 et la référence 0.025 sont
comparées aussi **entre solutions numériques**, afin de ne
pas confondre plancher spatial et convergence temporelle d'Euler explicite.
Les contrôles automatiques exigent une diminution des trois normes spatiales
et des différences RMS temporelles à la référence fine, mais aucun ordre fixé.
Les bilans indépendants et l'invariance axiale ont une tolérance de 1e-12,
conservatrice par rapport aux défauts mesurés de l'ordre de 1e-16.

## Résultats mesurés (Release, paramètres ci-dessus)

### Sensibilité à Transmissive en sortie subsonique

Sur R=1, les RMS (rho, rho*ur, rhoE) sont :

| Nr | rho | rho*ur | rhoE |
|---:|---:|---:|---:|
| 20 | 1.636569e-2 | 8.490812e-3 | 4.801464e-2 |
| 40 | 1.061558e-2 | 4.430478e-3 | 2.971967e-2 |
| 80 | 7.242448e-3 | 2.462733e-3 | 1.949779e-2 |

Le défaut de flux massique initial vaut -2.4786 %, -1.2447 %, -0.62369 %,
conforme à la projection du dernier anneau. Les ordres RMS finaux ne sont
que 0.552, 0.847 et 0.608. Ce résultat n'est pas présenté comme une
vérification d'ordre 1 indépendante de la fermeture de frontière.

Le contrôle étend le domaine jusqu'à R=2, sans modifier dx, dr ni l'état dans
la partie commune r<=1. Le solve étendu fournit les temps acceptés, auxquels
le solve initial est limité ; un sous-pas différent provoque un échec explicite.
Il n'y a ni solution analytique injectée dans les cellules, ni BC nouvelle.

| Nr | différence RMS rho | erreur rho R=1 | erreur rho domaine étendu, r<=1 |
|---:|---:|---:|---:|
| 40 | 1.054563e-2 | 1.063122e-2 | 2.571674e-3 |
| 80 | 7.232346e-3 | 7.249987e-3 | 1.295680e-3 |

Les mêmes diagnostics sont imprimés pour momentum et énergie. Il s'agit d'une
sensibilité au placement du bord, pas d'une erreur de BC isolée exactement.
En particulier, l'erreur de momentum sur R=1 bénéficie partiellement d'une
compensation : le domaine étendu ne réduit pas toutes les erreurs composante
par composante.

### Domaine complet avec sortie supersonique : R=2

Ce second domaine reste fixe pour ses trois raffinements. Le Mach analytique
au bord vaut `alpha0*R/a0 * s^(gamma-2)` ; il reste supérieur à 1.51 jusqu'à
t=0.2. Toutes les caractéristiques normales au bord extérieur sont donc
sortantes. Transmissive y est une approximation de trace sortante, sans
donnée entrante manquante. Ce choix évite d'introduire un mécanisme de BC
analytique dans la production ; il ne résout pas la limitation subsonique R=1.

| Nr | variable | L1 | RMS | Linf |
|---:|---|---:|---:|---:|
| 20 | rho | 5.936769e-3 | 1.084173e-2 | 7.468874e-2 |
| 40 | rho | 2.551438e-3 | 5.302104e-3 | 5.198699e-2 |
| 80 | rho | 1.051919e-3 | 2.436560e-3 | 3.049131e-2 |
| 20 | rho*ur | 1.751770e-2 | 1.826984e-2 | 3.260613e-2 |
| 40 | rho*ur | 8.395741e-3 | 8.804209e-3 | 1.695385e-2 |
| 80 | rho*ur | 4.029895e-3 | 4.232015e-3 | 9.332609e-3 |
| 20 | rhoE | 3.669953e-2 | 4.715205e-2 | 2.151994e-1 |
| 40 | rhoE | 1.696374e-2 | 2.239511e-2 | 1.524881e-1 |
| 80 | rhoE | 7.775413e-3 | 1.028219e-2 | 9.015929e-2 |

Ordres 40->80 (L1/RMS/Linf) : rho 1.278/1.122/0.770,
rho*ur 1.059/1.057/0.861, rhoE 1.125/1.123/0.758.
Les 140/281/563 pas sont acceptés sans rejet.
Le maximum de l'erreur de densité (Linf) se situe dans la première couronne.
Les trois niveaux ne démontrent pas un ordre uniforme 1 en Linf à l'axe ;
aucune explication asymptotique non vérifiée ni aucun seuil d'ordre artificiel
n'est utilisé pour effacer cette limite.

À Nr=80, RMS de pression par rapport au point exact : 3.311466e-3 ; par rapport
à l'EOS des moyennes exactes : 3.309488e-3. Ces diagnostics ne se substituent
pas aux normes conservatives.

### Convergence temporelle sur Nx=20, Nr=80, R=2

| CFL | pas | RMS différence à CFL=0.025 (rho / rho*ur / rhoE) |
|---:|---:|---|
| 0.4 | 141 | 3.806076e-4 / 9.315873e-4 / 2.506732e-3 |
| 0.2 | 282 | 1.775697e-4 / 4.345363e-4 / 1.169202e-3 |
| 0.1 | 563 | 7.610199e-5 / 1.862113e-4 / 5.010245e-4 |
| 0.025 | 2251 | référence |

Aucun rejet. Les ordres d'auto-convergence entre différences successives
0.4-0.2 et 0.2-0.1 valent 1.00073, 1.00116, 1.00127.
Les erreurs analytiques ne sont pas monotones en CFL : par exemple RMS(rhoE)
passe de 8.840881e-3 à 1.067072e-2 quand CFL descend de 0.4 à 0.025.
Il existe une compensation spatiale/temporelle ; une erreur analytique plus
petite à grand CFL n'indique donc pas un meilleur intégrateur.

## Bilans, tests et sorties

Les variations intégrales sont comparées aux flux réellement utilisés.
La source est reformée indépendamment depuis les états précédents avec
`G=2*pi*dx*dr` et l'EOS algébrique, pas depuis le diagnostic source du solveur.
Un test d'un pas mobile contrôle explicitement
`(Pr_new-Pr_old)/dt + B_radial = sum(G*p_EOS)` : omission ou inversion de
la source ne peut pas être masquée par un diagnostic modifié de la même façon.
Masse et énergie n'ont aucune source.

Les défauts sont normalisés par les intégrales initiales de masse et énergie,
et par `V*rho0*(alpha0*R+a0)` pour chaque momentum. Sur R=2/Nr=80, les défauts
globaux (M,Px,Pr,E) sont environ
`(-4.43e-18, 1.04e-19, 2.71e-17, -3.96e-17)`, et les défauts par pas <1.6e-16.
Les tests couvrent aussi repos, expansion faible/forte, échelles géométriques
0.01/1/100 (temps et alpha redimensionnés), raffinement et paramètres invalides.
Le calcul de référence ne change aucun composant de production.

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_VERIFICATION=ON
cmake --build build-release --target cfd_axisymmetric_euler_radial_expansion --parallel
./build-release/cfd_axisymmetric_euler_radial_expansion /tmp/radial-expansion
# Ouvrir /tmp/radial-expansion/radial_expansion.pvd dans ParaView.
```

Les quatre VTU (t=0,0.05,0.1,0.2) proviennent d'un solve séparé sur R=2/Nr=80,
dont les pas sont limités aux instants de sortie. Les tables de convergence
utilisent uniquement une limitation au temps final. Chaque sortie contient
rho, p (EOS numérique), ux, ur, vitesse, Mach, rhoE, les trois erreurs
conservatives et les pressions exactes ponctuelle/projetée, avec buffers réutilisés.

Limites restantes : frontière subsonique sans donnée entrante analytique ;
ordre Linf à l'axe non uniforme sur les niveaux étudiés ; précision/diffusion
du schéma Rusanov et d'Euler explicite d'ordre 1. Avant un cas radial subsonique
plus exigeant, une BC analytique de vérification ou une BC caractéristique
devra être étudiée séparément. Ni ces BC, ni la physique du lot suivant ne
sont développées ici. Aucune simulation RPA n'est lancée.
