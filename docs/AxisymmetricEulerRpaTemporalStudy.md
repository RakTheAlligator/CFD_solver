# RPA F3c — développement temporel contrôlé

Ce document conserve le compte rendu historique F3c : les mentions d'absence
de restart décrivent cette étape, pas le driver actuel. Le palier 3 a été
exécuté ultérieurement : 331720 pas sans rejet, 1 t_a atteint en 1186,08 s,
Mach col/sortie environ 1,0377/1,2181, toujours non stationnaire. Les fichiers
temporaires des anciens paliers ne sont plus disponibles ; ces chiffres sont
des repères rapportés, pas une nouvelle validation des fichiers disparus.
Le restart F3d et le stockage persistant sont décrits dans
[AxisymmetricEulerRestart.md](AxisymmetricEulerRestart.md).

| Palier | Cible (µs) | Atteint (µs) | Pas | Rejets | Mach col | Mach sortie | Mach max | Résidu max | Défaut conservatif normalisé max | Boucle (s) | Statut |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0,1 t_a | 50,925174 | 50,925174 | 27629 | 0 | 0,470780 | 0,780705 | 3,667353 | 8,91486e-3 | 3,18e-16 | 115,108 | TARGET_REACHED, non stationnaire |
| 0,5 t_a | 254,625869 | 254,625869 | 151155 | 0 | 0,824579 | 1,131875 | 5,713227 | 6,54301e-3 | 1,01e-15 | 737,108 | TARGET_REACHED, non stationnaire |
| 1,0 t_a (ultérieur) | 509,251737 | 509,251737 | 331720 | 0 | 1,0377 | 1,2181 | — | — | — | 1186,08 | TARGET_REACHED, non stationnaire |

Les Mach de section sont des estimations cellulaires pondérées par surface,
pas des traces exactes. Le maximum conservatif ci-dessus couvre les quatre
équations du ledger du driver ; le contrôle indépendant final donne
<=2,31e-16 au palier 1 et <=9,24e-16 au palier 2.
Atteindre le temps demandé n'est pas converger vers un état stationnaire.

## État et périmètre

Branche `feat/compressible-axisymmetric`, HEAD `f605229`. Le working tree
contenait déjà F3a et F3b non committés ; ces changements sont préservés.
F3c ne modifie que le driver RPA, ses outils/tests de vérification et ce
document. Aucun changement des flux, de l'EOS, du solveur Euler, de SIMPLE,
du maillage privé ou des BC. Pas de MUSCL, de Navier–Stokes ou de lot G.

Objectif : observer l'accélération et le démarrage du jet du modèle Euler
existant pendant un temps contrôlé, sans imposer l'étranglement ni chercher
un statut Converged artificiel.

## Configuration maintenue

- Maillage réel privé `cases/rpaNozzleSimple` : 2240 QUAD, 4675 faces,
  huit groupes, un domaine connecté. Le meshDict, le DXF et les autres fichiers
  privés sont contrôlés par empreintes avant/après ; aucune sortie dans ce cas.
- Gaz équivalent figé gamma=1,177, R=414,2 J/(kg K).
- Entrée : p statique=1170000 Pa, T=3343,5658 K, M=0,1508 ;
  TotalConditionsInlet avec P0=1185747,2501228245 Pa,
  T0=3350,294882116409 K, direction +x.
- Axis sur l'axe, SlipWall sur inner_wall/outer_wall/nozzle_lip ;
  Farfield existant sur upstream/radial/downstream.
- Pression ambiante nominale 1e-4 Pa ; T=300 K, ux=ur=0 et le même
  gaz équivalent extérieur restent des **hypothèses exploratoires**.
- Engine chaud uniforme au départ, quatre autres régions ambiantes selon
  RpaRegionLayout. 984 cellules Engine ; Band/Plume/Transition/Farfield
  contiennent respectivement 176/572/448/60 cellules.
- Rusanov interne, Euler explicite du premier ordre, CFL=0,01, sans plancher.
  EulerSteadyMonitor inchangé : résidus/variations/variation de débit 1e-7,
  déséquilibre de débit 1e-6, trois contrôles consécutifs après 100 pas.
  Observations toutes les 250 étapes et aux exports physiques ; aucun seuil
  de convergence n'est assoupli.

## Temps caractéristiques et budget annoncé

La définition F3b est conservée : longueur axiale de tuyère
L=0,650175702 m, son du gaz chaud initial
a=sqrt(gamma R T)=1276,7275086 m/s, donc
**t_a=L/a=5,09251737436e-4 s**. Le temps convectif basé sur la vitesse
d'entrée 192,5305083 m/s vaut 3,377000911e-3 s.

Ce n'est pas le temps acoustique de tout le domaine. À titre de références
cinématiques seulement, a_ambiant=382,4317194 m/s donne environ :

- 7,474 ms pour la distance sortie→aval, 2,85816514 m ;
- 3,363 ms pour l'écart radial sortie→farfield supérieur, 1,286174313 m ;
- 1,700 ms pour la distance axiale sortie→frontière amont, 0,650175702 m.

Ces divisions par une vitesse figée ne bornent pas l'arrivée des ondes non
linéaires dans un milieu qui accélère. Le temps d'établissement stationnaire
est inconnu et ne peut être déduit d'un seul temps de traversée de la tuyère.

À partir de F3b : dt moyen=1,93547e-9 s, coût de boucle=3,781 ms/pas.
L'estimation initiale pour le premier palier était 26312 pas / 99,5 s.
**Budget explicite appliqué : 100000 pas, 300 s de boucle, 20 rejets cumulés.**
La limite murale est vérifiée entre mises à jour, sans interruption au milieu
d'un flux/export ; ce n'est pas une limite temps réel stricte du processus.

Mesuré au premier palier :

- 27629 pas, dt moyen=1,8431783e-9 s, maximum=1,9468741e-9 s ;
- minimum=7,2237165e-10 s, dû à un export physique ; en excluant les cinq
  pas tronqués aux exports : minimum=1,7670022e-9 s ;
- CFL configuré constant 0,01, reporté comme borne supérieure, pas comme
  mesure indépendante du CFL effectif après troncature ; aucun rejet ;
- 115,108 s de boucle, 116,66 s murales du processus, CPU utilisateur
  116,94 s et système 0,41 s (incluant les sous-processus) ;
- RSS maximale rapportée par `/usr/bin/time -v` : 51604 KiB ;
- sorties CFD brutes (VTU/PVD/CSV/métadonnées, avant rendus) : 16487467 octets,
  dont six VTU d'environ 0,91 Mo chacun.

Les coûts observés par quart sont 4,532 / 4,234 / 3,675 / 4,228 ms/pas.
Ils ne sont pas constants. Des validations s'exécutaient aussi sur la machine :
il ne s'agit pas d'un benchmark CPU isolé et aucune variation de coût n'est
attribuée exclusivement à la physique ou au Farfield.

## Exécution et continuation

Le driver ne dispose pas d'un restart conservatif. Un VTU n'en est pas un :
temps, compteurs, ledger cumulé et état du moniteur ne seraient pas restaurés
par le simple rechargement de ses champs. Le premier palier repart donc de t=0.
Aucune API de restauration ni modification de production n'est ajoutée.

L'option finale `--physical-exports` conserve les appels F3b existants sans
changement. Elle remplace les exports par pas par cinq cibles physiques
(0,1 ; 0,25 ; 0,5 ; 0,75 ; 1) fois MAX_TIME, plus l'état initial.
Pour le palier 1 cela donne t/t_a=0 ; 0,01 ; 0,025 ; 0,05 ; 0,075 ; 0,1.
`solver.step(next_target)` tronque le pas ; tous les temps écrits sont ceux
réellement acceptés, sans réétiquetage ni interpolation temporelle des champs.
Cela introduit cinq pas tronqués, pas une modification du CFL demandé.

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_VERIFICATION=ON -DBUILD_TESTING=ON \
  -DPython3_EXECUTABLE="$HOME/.venvs/cfd-rpa/bin/python"
cmake --build build-release --target cfd_axisymmetric_euler_rpa_exploratory

# Choisir une sortie neuve : un dossier non vide est refusé, jamais nettoyé.
/usr/bin/time -v -o /tmp/cfd-f3c-stage1-resources.log \
  ./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/cfd-rpa-f3c-stage1 \
  300 0 0.01 100000 10000 0 0.0000509251737436 250 20 300 \
  --physical-exports
python3 verification/validate_rpa_euler_output.py /tmp/cfd-rpa-f3c-stage1
```

Les arguments restent CASE OUTPUT T_inf ux_inf CFL MAX_STEPS SAVE_EVERY
ur_inf MAX_TIME DIAG_EVERY MAX_REJECTED MAX_WALL_SECONDS, suivis du flag
facultatif, puis éventuellement un unique `COMPARISON_TIME` strictement
intérieur à l'intervalle. Ce temps est inséré dans l'ordre croissant sans
dupliquer une cible existante ; les cinq cibles par défaut restent inchangées.
SAVE_EVERY est ignoré en mode physique. `TARGET_REACHED` signifie
que MAX_TIME a été atteint, pas que le moniteur est convergé. STEP_LIMIT,
WALL_LIMIT et REJECT_LIMIT restent des paliers incomplets si la cible manque.
`status.txt` fournit aussi `target_reached` et la cause/phase des exceptions.
Les diagnostics et le dernier état accepté restent disponibles en cas d'échec.

## Sections réelles et diagnostics

Le nouvel outil de vérification sélectionne **une fois** les faces réellement
alignées sur x_inlet, x_throat et x_exit issus du domaine CAD discrétisé.
Chaque section comporte huit anneaux, avec couverture radiale sans trou ni
recouvrement. Une section non alignée ou partielle est rejetée ; aucun choix
silencieux de « cellule la plus proche ». Les aires proviennent des métriques
de révolution existantes, pas des longueurs méridiennes.

Les scalaires Mach total, Mach axial ux/a, p et T sont des moyennes pondérées
par aire des cellules Engine adjacentes : deux côtés au col, côté Engine
seul à l'entrée/sortie. Les min/max décrivent ces estimations par anneau,
pas des traces Riemann. Les distributions sont également accessibles dans
les champs cellulaires exportés.

Les **débits** sont au contraire des sommes de flux numériques conservatifs
intégrés vers +x : Rusanov exact du code aux faces internes, flux caractéristique
existant à l'entrée. Ils ne sont pas le produit de moyennes de rho et de vitesse.
Il n'y a aucune allocation par évaluation de section et aucun recalcul de sa
connectivité pendant les observations.

`sections.csv` contient ces valeurs. `cell_diagnostics.csv` ajoute t/t_a, dt,
bornes rho/p, min/max T, Mach max global/Engine, débits, résidu max et temps mural.
Pour le palier 2, la masse Engine est aussi accumulée pendant le même parcours
cellulaire et écrite dans `engine_mass`, sans nouvelle traversée du maillage.
`convergence.csv` conserve les quatre résidus, variations et critères existants.
`steps.csv` conserve à chaque pas accepté les intégrales, flux et défauts de
bilan. `boundary_groups.csv` ajoute min/max des pressions des traces à ses
régimes, Mach normaux et débits courants. Ces derniers ne sont pas les flux
de l'étape précédente utilisés dans le ledger temporel.

## Résultats physiques calculés

| t/t_a | Pas | Mach col | Mach sortie | p col (Pa) | p sortie (Pa) | T col (K) | T sortie (K) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0,150800 | 0,150800 | 1170000 | 1170000 | 3343,566 | 3343,566 |
| 0,01 | 2631 | 0,172810 | 0,302153 | 1161571 | 920827 | 3340,503 | 3329,851 |
| 0,025 | 6623 | 0,237389 | 0,472756 | 1155043 | 734054 | 3339,258 | 3285,308 |
| 0,05 | 13458 | 0,343028 | 0,629347 | 1119956 | 607611 | 3327,509 | 3211,958 |
| 0,075 | 20482 | 0,413333 | 0,717962 | 1050610 | 544957 | 3300,926 | 3152,377 |
| 0,1 | 27629 | 0,470780 | 0,780705 | 995651 | 501875 | 3280,044 | 3102,079 |

À la fin, Mach axial moyen=0,470071 au col et 0,774735 à la sortie ;
les estimations par anneau varient de 0,460011 à 0,487094 au col et de
0,747552 à 0,836918 à la sortie. Mach max Engine=0,836918 : **aucun
étranglement ni région supersonique interne démontré** à ce palier.
Mach max global=3,667353 dans TRANSITION vers (x,r)=(1,15728 ; 0,41903) m.
Les régions géométriques ne sont pas des espèces ou des fractions de mélange.

| Région finale | Masse (kg) | Mach max | Pression max (Pa) |
|---|---:|---:|---:|
| Engine | 8,786421e-3 | 0,836918 | 1290624 |
| Band | 1,298930e-5 | 3,361843 | 8945,109 |
| Plume | 1,833938e-3 | 3,609765 | 431979,665 |
| Transition | 6,283178e-5 | 3,667353 | 6393,843 |
| Farfield | 1,160462e-8 | 0,218479 | 1,495793e-4 |

La région volumique Farfield n'est pas le groupe de faces extérieures : son
Mach cellulaire maximal ne classe pas le régime des traces de frontière.

Débits instantanés vers +x : entrée 0,417436 kg/s, col 0,290864 kg/s,
sortie Engine→Plume 30,772445 kg/s. Ce dernier ne constitue surtout pas un
débit moteur stationnaire : il était déjà 45,037696 kg/s dans le flux initial
chaud/ambiant discontinu, atteint 49,049552 au premier export, puis décroît.
Le stock Engine passe de 0,010674916 à 0,008786421 kg (−17,69 %).
L'écoulement évacue donc une partie du gaz arbitrairement chaud initialement
réparti dans le divergent ; l'injection amont seule n'alimente pas ce débit.

Le calcul reste positif à chaque état accepté : p_min=1e-4 Pa,
rho_min=8,04764204e-10 kg/m3. Maxima finaux p=1,290624 MPa,
rho=0,918146 kg/m3 ; maxima sur l'historique accepté p=1,362854 MPa,
rho=0,961529 kg/m3. Températures observées entre 300 K (à l'arrondi près)
et 3421,959 K ; maximum final 3393,734 K. Ce démarrage non stationnaire
comprend aussi des compressions dues à l'initialisation non tangente aux parois.

Les vues montrent l'extension de la détente et du domaine extérieur
supersonique. Comme marqueurs **opérationnels**, les centres extérieurs avec
p>2 p_inf / 10 p_inf / 100 p_inf atteignent respectivement x=1,50965 /
1,50952 / 1,39895 m et r=0,81270 / 0,81057 / 0,70452 m. Ces contours
dépendent du seuil et de la diffusion du schéma : ce ne sont ni des fronts
physiques exacts ni une démonstration de choc. Aucun test Rankine–Hugoniot
ni étude de résolution n'a été effectué pour identifier une structure de choc.

## Farfield, résidus et conservation

| Groupe | Débit sortant signé (kg/s) | Sub entrant/sortant | Sup entrant/sortant | Mach normal max | p trace min/max (Pa) |
|---|---:|---:|---:|---:|---:|
| upstream | 2,43463e-12 | 1 / 10 | 0 / 0 | 1,05561e-5 | 1e-4 / 1,00001242468e-4 |
| radial | 2,14306e-12 | 1 / 19 | 0 / 0 | 2,72880e-6 | 1e-4 / 1,00000321180e-4 |
| downstream | 3,20768e-24 | 5 / 11 | 0 / 0 | 2,21974e-17 | 1e-4 / 1e-4 |

La classification est celle de la trace normale Farfield existante, pas du
Mach total. L'entrée reste huit faces subsoniques entrantes à M_n=−0,1508.
Les minima M_n des trois farfields sont respectivement −1,13e-17,
−1,18e-18 et −1,40e-17. Ces entrées minuscules ne prouvent pas un entraînement
physique ; en revanche les petites sorties amont/radiales dépassent désormais
le seul bruit d'arrondi F3b. Elles signalent une faible perturbation numérique,
pas une plume établie ou une arrivée de front physique identifiée. Le débit
aval reste au niveau de l'arrondi. Débit extérieur total=4,57769e-12 kg/s.

Résidus finaux normalisés (rho, rhoux, rhour, rhoE) :
(1,675213e-3 ; 1,467484e-3 ; 8,366936e-4 ; 8,914861e-3).
Depuis F3b à 3000 pas, trois diminuent mais le radial augmente ; aucune
convergence simultanée. Déséquilibre massique entrée/sortie=0,999999999978.
La masse totale augmente conformément à l'injection moins les sorties :
un tel déséquilibre transitoire n'est pas une erreur de conservation.
Les grandes échelles chaudes de F1 sur tout le domaine restent peu sensibles
aux erreurs relatives locales de l'ambiance très diluée.

Le ledger du driver vérifie les flux sur toutes les frontières, l'accumulation
et la source radiale. Ses défauts absolus finaux sont
(5,66e-18 kg ; −1,70e-15 kg m/s ; 8,93e-17 kg m/s ; 5,32e-11 J).
Les défauts relatifs à la somme des intégrales et flux intégrés sont tous
<=3,18e-16. Le contrôle indépendant reconstruit les volumes depuis les
polygones VTU, somme avec `math.fsum` et intègre l'historique des taux ;
il obtient <=2,31e-16. Les différences absolues restantes sont compatibles
avec l'arrondi et les soustractions des intégrales. Les bornes restent
proportionnelles aux échelles effectivement sommées et au nombre d'opérations,
pas à un seuil physique large choisi après le calcul.

## Temps et sorties visuelles

| Compteur | Secondes |
|---|---:|
| Assemblage production | 97,3281 |
| Toutes BC, sous-ensemble de l'assemblage | 42,2593 |
| Mise à jour conservative | 16,4199 |
| Monitoring, incluant des préparations d'assemblage | 0,62523 |
| Farfield échantillonné, inclus dans monitoring | 0,16309 |
| Sections géométriques et CSV, inclus dans monitoring | 0,00715 |
| Export, état initial inclus | 0,17139 |
| Boucle | 115,1085 |

Ne pas additionner ces compteurs recouvrants. Le temps BC de production
regroupe aussi entrée et parois ; il n'isole pas le Farfield à chaque pas.
Le timer Farfield diagnostique ne mesure que les observations échantillonnées.
Aucun nouveau timer ni optimisation de production n'est introduit.

PVD : `/tmp/cfd-rpa-f3c-stage1/rpa.pvd`. Les six VTU conservent rho, p,
T, ux, ur, vitesse, Mach, rhoE, log10_p_Pa, log10_rho, rpa_region et U.
Le validateur indépendant contrôle XML, connexions/orientation, volumes,
positivité, EOS, états initiaux, chronologie, temps acceptés, cibles physiques,
bilans et moyennes de section recalculées depuis les VTU.
Les débits numériques internes col/sortie sont également recalculés
indépendamment à partir des primitives exportées et du flux scalaire Rusanov.
PVDReader a chargé les six instants. Les rendus utilisent des plages fixes :
Mach [0,4], log10(p/Pa) [−4,7], T [300,3500] K, vitesse [0,4200] m/s.
Les PNG sont annotés avec le vrai temps et « NOT STEADY » ; l'animation
`Mach_temporal.mp4` conserve des durées proportionnelles aux intervalles
physiques, avec un arrêt final d'une seconde. Aucun affichage GUI interactif
n'est revendiqué pour ce palier ; les rendus ParaView ont été examinés.

```sh
paraview /tmp/cfd-rpa-f3c-stage1/rpa.pvd
```

L'environnement lancé depuis Snap peut demander un nettoyage **local au
processus** des variables GTK_PATH/GTK_MODULES/GIO_MODULE_DIR ; aucune
bibliothèque système n'a été changée. Tous les résultats F3b sont préservés.

## Validation logicielle et décision

Les tests courts ajoutés couvrent le temps acoustique, les cibles invalides
ou non représentables, leur ordre, la troncature à un état réellement accepté,
la limite de pas, les sections complètes/partielles, la moyenne côté Engine
et un débit analytique uniforme avec le signe d'entrée correct. Aucun cas
privé ni calcul long n'est ajouté à CTest. Des tests CLI réels séparés ont
confirmé TARGET_REACHED (5 pas très courts), STEP_LIMIT (1 pas, cible manquée)
et WALL_LIMIT (0 pas, cible manquée), avec sorties indépendamment contrôlées.

Builds complets Debug/Release et ASan+UBSan réussis ; CTest complet **54/54**
dans chacune des trois configurations. Le nouveau mode physique a également
passé un court run privé instrumenté (8 pas, six exports), puis son contrôle
indépendant. LeakSanitizer reste désactivé (`detect_leaks=0`), conformément
à la limitation Gmsh déjà documentée : aucune vérification des fuites revendiquée.
clang-format et clang-tidy strict ciblés passent ; un accès optional non
reconnu comme gardé par clang-tidy a été rendu explicitement gardé, sans
suppression de diagnostic. `git diff --check` passe. Les empreintes des neuf
fichiers privés, résultats antérieurs compris, restent identiques.

**Prolongation possible sous conditions**, sans modification physique :
positivité, flux Farfield calculables et bilans sont satisfaisants. Le col
accélère mais reste subsonique ; le démarrage est dominé par le stock initial
et ne valide pas encore le fonctionnement nominal RPA.

Projection initiale après le palier 1, pas résultat : à dt moyen observé,
le palier 2 demanderait environ
138145 pas et 576 s depuis zéro ; au dt de fin, environ 144000 pas. Prévoir
plutôt 10–15 min et une limite explicite de 900 s / 500000 pas, **à confirmer
par l'utilisateur avant lancement**. Le coût peut augmenter ; le palier 3
serait au moins de l'ordre de 20 min et exige une nouvelle décision après
analyse du second. L'autorisation a ensuite été donnée pour le palier 2,
avec ces plafonds ; ses résultats réels sont détaillés ci-dessous.
Le palier 3 reste non exécuté et nécessite un nouvel accord.

Un nouveau lancement sans restart recalculerait les 115 s du premier palier
(environ 20 % de l'estimation totale du second). Une reprise exacte impose
de traiter temps/compteurs/ledger/moniteur, pas seulement u/rho/E d'un VTU.
Ce coût évitable est identifié, mais ne justifie pas à lui seul d'introduire
dans ce lot une nouvelle API de restauration de production non validée.

La validité d'Euler dans l'ambiance raréfiée reste non démontrée. Ni les
bilans, ni les contours de Mach, ni ce palier stable ne valident une plume
réelle à 200–800 km. Pas de chimie, espèces, viscosité, turbulence ou modèle
cinétique ; aucune performance spatiale n'est revendiquée. Aucun commit/push.

## Palier 2 autorisé : résultat à 0,5 t_a

### Protocole et modification limitée des diagnostics

Calcul entièrement repris de t=0, sans restart, à configuration physique,
maillage et CFL=0,01 inchangés. Limites autorisées : 900 s de boucle,
500000 pas acceptés, même plafond de 20 rejets cumulés. Le calcul a atteint
la cible avant ces plafonds : **151155 pas, 0 rejet, 737,108 s de boucle**,
statut **TARGET_REACHED**. Le statut du moniteur est TimeLimit, pas Converged.

Le seul ajout au plan d'exports est un temps de comparaison facultatif ;
ici c'est le temps réellement accepté au palier 1, 5,09251737436e-5 s.
Les sept états exportés sont à t/t_a=0 ; 0,05 ; 0,1 ; 0,125 ; 0,25 ;
0,375 ; 0,5. La masse Engine est ajoutée au diagnostic déjà parcouru,
sans allocation ni parcours de cellules supplémentaire par observation.
Aucune méthode numérique ou option du solveur de production n'a changé.

```sh
/usr/bin/time -v -o /tmp/cfd-f3c-stage2-resources.log \
  ./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/cfd-rpa-f3c-stage2 \
  300 0 0.01 500000 5000 0 0.000254625868718 250 20 900 \
  --physical-exports 0.0000509251737436
python3 verification/validate_rpa_euler_output.py /tmp/cfd-rpa-f3c-stage2
```

Pour reproduire, choisir une sortie neuve ; les dossiers non vides restent
refusés. Les empreintes des **41 fichiers du palier 1**, ainsi que des neuf
fichiers du cas privé, restent inchangées. Les nouvelles vues comparatives
sont écrites séparément, jamais dans le dossier du palier 1.

Au temps commun 0,1 t_a, les différences max entre les deux lancements sont
3,17e-10 kg/m3 en rho, 0,0005001 Pa en p, 0,000156 m/s en ux,
0,000102 m/s en ur et 1,36e-7 en Mach. L'écart de Mach moyen au col
est 6,71e-11. Les instants correspondent exactement, mais les troncatures
antérieures aux exports sont différentes : l'identité bit à bit n'est donc
ni attendue ni revendiquée. La trajectoire calculée est cohérente à cette
précision sans changement de formulation.

### Débits et vidange de l'inventaire initial

| t/t_a | Mach col | Mach sortie | Q entrée (kg/s) | Q col (kg/s) | Q sortie (kg/s) | Masse Engine (g) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0,150800 | 0,150800 | 0,417436 | 0,104359 | 45,037696 | 10,674916 |
| 0,05 | 0,343028 | 0,629347 | 0,417436 | 0,236966 | 34,671570 | 9,601422 |
| 0,1 | 0,470780 | 0,780705 | 0,417436 | 0,290864 | 30,772445 | 8,786421 |
| 0,125 | 0,522082 | 0,830102 | 0,417436 | 0,312436 | 29,603885 | 8,407642 |
| 0,25 | 0,664806 | 0,988852 | 0,417397 | 0,365017 | 25,120397 | 6,697350 |
| 0,375 | 0,708118 | 1,076140 | 0,341686 | 0,377570 | 20,579360 | 5,265548 |
| 0,5 | 0,824579 | 1,131875 | 0,121987 | 0,402073 | 16,307050 | 4,109974 |

La masse Engine baisse de **61,50 %** depuis l'initialisation. Le flux de
sortie diminue de **47,01 %** par rapport au palier 1, mais reste 39,1 fois
le débit initialement injecté et 133,7 fois le débit entrant courant.
La vidange du gaz initialement chaud dans le divergent reste dominante ;
ce flux n'est ni un débit moteur stationnaire ni une prédiction de performance.

**0,4174 kg/s n'est pas une consigne de débit.** P0 et T0 restent imposées,
mais TotalConditionsInlet détermine le flux avec l'état intérieur. Le débit
entrant finit à 0,121987 kg/s, après un minimum observé de 0,111048 kg/s
et un léger rebond. La trace d'entrée reste subsonique entrante, à Mach
normal environ −0,0435715 ; sa pression statique est 1184423,32 Pa.
Cette évolution et ce rebond témoignent d'un démarrage encore fortement
instationnaire, pas d'une modification des conditions nominales imposées.

### Col, sortie et développement extérieur

À 0,5 t_a, les estimations par anneau au col vont de Mach 0,786394
à 0,871108 ; les Mach axiaux vont de 0,786238 à 0,868204. Les estimations
au col restent donc toutes subsoniques : **aucun étranglement stationnaire
n'est démontré**. À la sortie, les valeurs Engine adjacentes sont supersoniques
sur les huit anneaux : Mach total 1,086057–1,195681, axial 1,083341–1,172687.
Le Mach global 5,713227 est situé dans PLUME près de
(x,r)=(2,89869 ; 0,0153792) m, pas au col.

Pressions moyennes col/sortie : 744363,65 / 196187,70 Pa ; températures :
3185,827 / 2622,350 K. Malgré l'accélération, le divergent n'a pas atteint
les conditions stationnaires RPA. Les comparaisons ponctuelles à RPA ne
constitueraient pas ici une vérification spatiale ou thermochimique.

Le jet numérique s'étend nettement dans le domaine extérieur. Les masses
finales des régions sont respectivement Engine 4,109974e-3 kg,
Band 1,340194e-5 kg, Plume 3,810153e-3 kg, Transition 2,560817e-3 kg,
Farfield 2,440446e-4 kg. Les maxima de Mach correspondants sont
1,195681 / 3,432308 / 5,713227 / 5,500538 / 5,156274.
Il s'agit d'étiquettes géométriques, pas d'espèces. Les structures visibles
n'ont pas été identifiées comme des chocs à l'aide de Rankine–Hugoniot.

### Farfield effectivement perturbé

| Groupe à 0,5 t_a | Débit sortant (kg/s) | Sub entrant/sortant | Sup entrant/sortant | Mach normal min/max | Pression trace min/max (Pa) |
|---|---:|---:|---:|---:|---:|
| upstream | 0,0596610 | 0 / 0 | 0 / 11 | 1,41678 / 3,43225 | 1,12160 / 6,00773 |
| radial | 0,7190147 | 0 / 2 | 0 / 18 | 0,060078 / 3,95756 | 1,07637e-4 / 18,49577 |
| downstream | 2,59100e-6 | 0 / 7 | 0 / 9 | 0,074847 / 2,01110 | 1,09608e-4 / 1,81084e-3 |

Les trois groupes sont désormais perturbés, amont et radial bien davantage
que l'aval en débit. Le débit extérieur total atteint 0,778678 kg/s,
contre une entrée de 0,121987 kg/s : à l'instant final la masse globale
diminue, avec dM/dt=−0,656692 kg/s. Le déséquilibre de débit 0,729119 est
compatible avec ce transitoire, sans être un défaut de conservation.

Les pressions des traces sortantes ne sont pas tenues égales à p_inf :
Farfield n'est pas PressureOutlet. L'état extérieur fourni garde bien sa
pression nominale 1e-4 Pa ; le flux et sa classification résultent du même
problème normal de Riemann que dans l'assemblage. Aucune entrée supersonique
n'est observée à l'instant final. Ces perturbations numériques ne démontrent
pas à elles seules une position de front physique précise dans le quasi-vide.

### Positivité, résidus et bilans

Sur les 151155 états acceptés : p_min=1e-4 Pa,
rho_min=8,0446872e-10 kg/m3, p_max=1,3628543 MPa,
rho_max=0,9615292 kg/m3. Les températures échantillonnées restent proches
de 300 K au minimum et atteignent 3421,976 K au maximum. Les minima finaux
sont p=1,1144745e-4 Pa, rho=8,1891157e-10 kg/m3, T=328,5663 K.
La positivité est préservée sans plancher ni rejet de pas.

| t/t_a | R rho | R rhoux | R rhour | R rhoE | Défaut conservatif normalisé max indépendant |
|---:|---:|---:|---:|---:|---:|
| 0 | 2,41784e-3 | 1,73294e-3 | 2,63789e-4 | 1,18755e-2 | 0 |
| 0,05 | 1,85827e-3 | 1,58612e-3 | 6,21217e-4 | 9,74437e-3 | 4,33e-16 |
| 0,1 | 1,67521e-3 | 1,46748e-3 | 8,36694e-4 | 8,91486e-3 | 6,53e-16 |
| 0,125 | 1,63575e-3 | 1,42132e-3 | 9,19238e-4 | 8,67718e-3 | 9,59e-16 |
| 0,25 | 1,50036e-3 | 1,29958e-3 | 1,09342e-3 | 7,97306e-3 | 9,54e-16 |
| 0,375 | 1,38259e-3 | 1,40897e-3 | 1,21356e-3 | 7,27395e-3 | 9,27e-16 |
| 0,5 | 1,26358e-3 | 1,34416e-3 | 1,22770e-3 | 6,54301e-3 | 9,23e-16 |

Les résidus massique et énergétique diminuent mais le radial augmente ;
le momentum axial n'est pas monotone. Aucun critère simultané stationnaire
n'est satisfait. Les valeurs absolues finales sont
(47,1833 kg/s ; 64081,56 N ; 58529,69 N ; 3,98253e8 W).

Défauts absolus cumulés du driver : masse 1,03e-17 kg, momentum axial
1,17e-14 kg m/s, radial 3,21e-15 kg m/s, énergie 1,70e-10 J.
Les défauts normalisés du driver sont <=1,01e-15. Le contrôle indépendant
depuis les polygones et l'historique retrouve <=9,24e-16 au dernier état
avec une normalisation incluant les taux effectivement sommés. La source
radiale est incluse. Les différences d'arrondi entre les deux calculs
ne sont pas des pertes physiques de conservation.

### Coût réel et limites

dt moyen=1,684535e-9 s. Minimum hors troncatures d'export=1,550351e-9 s,
maximum=1,946874e-9 s. Minimum absolu=1,129684e-10 s, à un export tronqué.
La borne CFL reste exactement 0,01 ; les six troncatures ne sont pas une
modification de cette option. Coût moyen de boucle=4,877 ms/pas ; par
quart approximatif : 5,829 / 4,939 / 4,924 / 3,811 ms/pas. Les validations
parallèles du début affectent ces temps, qui ne constituent pas un benchmark
CPU isolé. Aucune hausse ou baisse n'est attribuée sans preuve au seul Farfield.

Temps mur du processus 739,06 s, CPU utilisateur/système 731,85 / 4,70 s,
RSS max 51604 KiB. Sorties CFD brutes avant rendus : 68217669 octets.
Assemblage=622,780 s, toutes BC incluses=251,114 s, mise à jour=100,562 s,
monitoring=4,0164 s, Farfield diagnostique inclus=0,9559 s,
sections incluses=0,0520 s, export=0,2669 s. Comme au premier palier,
les compteurs recouvrants ne doivent pas être additionnés.

### Contrôles et livrables

L'ancien validateur a signalé à tort un incrément temporel au pas 144441 :
soustraire les temps proches donnait une erreur de 2,565e-20 s, inférieure
à une demi-ULP du temps courant (ULP=5,421e-20 s), mais supérieure à son
ancien seuil relatif à dt. Les 151155 lignes satisfont toutes exactement
`time == previous_time + dt` après relecture des 17 chiffres CSV.
Le contrôle compare maintenant cette accumulation et exige une progression
strictement positive. Aucun seuil du solveur ou du résidu n'est relâché ;
un test reproduit l'annulation et vérifie le rejet de timestamps corrompus
d'une ULP. Les fichiers calculés n'ont pas été modifiés pour passer ce contrôle.

Le validateur complet passe les sept snapshots, leur EOS, leurs volumes,
les deux débits numériques internes, la masse Engine, les cibles et les quatre
bilans. Les six anciens snapshots du palier 1 passent aussi le nouveau contrôle.
CTest Debug, Release et ASan+UBSan : **54/54** dans chaque configuration ;
clang-format, clang-tidy strict ciblé et `git diff --check` passent.
Les tests incluent le tri, la déduplication et le rejet des temps de
comparaison invalides. LeakSanitizer reste désactivé, sans revendication de
contrôle des fuites.

PVD : `/tmp/cfd-rpa-f3c-stage2/rpa.pvd`. PVDReader charge les sept instants.
Les rendus des **deux paliers** sont dans `/tmp/cfd-rpa-f3c-stage2-views`,
avec mêmes caméras et plages fixes partagées : Mach [0,6], log10(p) [−4,7],
T [200,3500] K, vitesse [0,5600] m/s. Ces plages sont seulement visuelles.
L'animation `Mach_temporal.mp4` est annotée avec les temps réellement acceptés ;
aucun résultat du palier 1 n'est réécrit. Les distributions radiales
recalculées indépendamment sont disponibles dans
`/tmp/cfd-rpa-f3c-stage2-analysis/section_radial_profiles.csv`.
Une seconde animation `Mach_domain.mp4` montre le domaine entier, avec la
même plage Mach [0,6], afin de ne pas masquer les interactions avec les
frontières dans les vues rapprochées.

Le palier 2 est donc **numériquement admissible et conservatif jusqu'à sa
cible**, mais reste dominé par l'inventaire initial et un débit d'entrée
non stabilisé. Les frontières extérieures participent maintenant au transitoire.
Il ne valide ni une plume réelle raréfiée ni les performances nominales RPA.
Le palier 3 n'a pas été lancé ; un nouvel accord et un examen de ces résultats
sont nécessaires avant toute prolongation. Aucun commit, aucun push.
