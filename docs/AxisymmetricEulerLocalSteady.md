# Euler axisymétrique stationnaire — pseudo-temps local explicite

Cette V1 recherche une racine du **résidu FV discret**. Elle ne modélise pas
l'évolution physique du jet et ne garantit pas l'existence d'un état stationnaire.
Flux Rusanov, gaz parfait, conditions aux limites et source radiale sont ceux du
solveur transitoire, sans approximation supplémentaire ni plancher de positivité.

## Interfaces et séparation des modes

| Propriété | Transitoire | Stationnaire |
| --- | --- | --- |
| Type | `AxisymmetricEulerSolver` | `AxisymmetricEulerSteadySolver` |
| Objectif | Évolution physique | Résoudre `R(U)=0` |
| Avancement | `step(final_time)` | `iterate()` / `solve()` |
| Pas | Global, limité par toutes les cellules | Local, à partir des mêmes sommes spectrales |
| Compteur | Pas acceptés et temps physique | Itérations acceptées, sans temps physique |
| Conservation | Bilan temporel avec flux intégrés | Résidus et bilans instantanés flux/source |
| Historique de reprise | Checkpoints transitoires existants | Non implémenté dans cette V1 |
| Observation | `EulerSteadyMonitor` du transitoire | Contrôles propres au solveur local |

Le nouveau type possède un noyau `AxisymmetricEulerSolver` privé. Une amitié
limitée permet de réutiliser son assembleur et ses buffers, sans appeler son
avancement physique ni exposer ses diagnostics/checkpoints temporels. Aucun enum
dans les options historiques, aucun changement de comportement par défaut.
`EulerSteadyMonitor` reste inchangé : sa dépendance au temps physique et sa
normalisation par le seul volume global ne conviennent pas directement au mode local.

Le `Mesh` doit rester vivant et immobile. Le champ initial est copié profondément
une fois. Les vues de l'état et des pas locaux sont non propriétaires et peuvent
changer à l'itération suivante. Les étiquettes régionales facultatives sont copiées
et remappées une fois ; elles n'imposent aucune topologie F4.

## Formulation et rejet transactionnel

Avec les orientations owner/neighbor existantes :

```
R_C = -sum_f(sign_C,f * S_f * Fhat_f) + (0,0,G_C*p_C,0)
G_C = 2*pi*A_2D,C
wave_C = sum_f(S_f * lambda_f)
dtau_C = CFL_steady * V_C / wave_C
Ucandidate_C = Uold_C + (alpha*dtau_C/V_C)*R_C(Uold)
```

`S_f` est la norme de la surface axisymétrique. `lambda_f` est la vitesse réellement
retournée par le flux interne ou la condition aux limites ; ce n'est pas un second
estimateur approximatif par cellule. Les faces d'axe de surface nulle sont sautées
par l'assembleur existant. Aucun calcul ne divise par leur surface.

Toutes les cellules lisent le même état précédent : mise à jour **Jacobi**, dans
le buffer candidat distinct. Densité, énergie interne, pression, température et
représentabilité sont validées par l'EOS existante. Si un candidat échoue,
`alpha=1,1/2,1/4,...` est réduit globalement et **tout** le candidat est recalculé
depuis le dernier état accepté. Aucune cellule active n'est modifiée avant
acceptation. Avant la permutation, chaque face frontière de surface non nulle
est évaluée sur le candidat avec le même dispatch de flux, les mêmes normales,
BC et aires de groupe que l'assemblage. Une `EulerCharacteristicCompatibilityError`
(branche/régime caractéristique incompatible) rejette ce candidat et utilise la
même réduction globale. Les autres erreurs d'évaluation des BC (données,
géométrie ou arithmétique) restent fatales ; l'assemblage de l'état courant n'est jamais
transformées en rejets. Un état initial incompatible échoue avant toute tentative.
Pas nul, spectral nul/non fini, overflow/underflow du pas
produisent un diagnostic explicite, sans epsilon ou plafonnement caché.

En arithmétique exacte, les facteurs `alpha*dtau_C/V_C` strictement positifs
préservent les racines `R=0` : seule la trajectoire itérative change. En flottant,
une mise à jour trop petite peut ne plus modifier le champ alors que `R` reste
grand ; c'est précisément pourquoi le résidu explicite reste obligatoire.

Le CFL des faces n'est pas un théorème de positivité pour la source axisymétrique.
Les tentatives et les rejets cumulés sont bornés. L'absence de rejet ne démontre
ni la précision ni la convergence stationnaire.

Les buffers de pas, états, résidus, comparaison des champs et sommes régionales
sont persistants. Une liste des faces frontières non nulles est construite une
fois : chaque tentative valide ajoute O(Nfrontières) évaluations, sans second
assemblage des faces internes ni nouvelle allocation sur le chemin sans rejet.
Les messages de rejet/échec peuvent allouer. La validation ne modifie ni l'état
accepté ni les caches ; ceux-ci sont invalidés seulement à l'acceptation.
Elle garantit la réévaluation des BC immuables sur le champ accepté, pas la
représentabilité de toutes les sommes du prochain assemblage, ni l'existence
d'une solution stationnaire. Une approche de borne à pas presque nul reste
soumise aux résidus explicites, limites de tentatives/rejets, stagnation et temps.
Aucun cache thermodynamique,
aucune optimisation du flux, aucun parallélisme n'est introduit.

## Résidus, échelles et convergence

Références positives, finies et fixes, configurables dans `EulerLocalSteadyOptions` :

```
Qref = (rho_ref, rho_ref*Uref, rho_ref*Uref, rho_ref*Uref^2)
tref = Lref/Uref
```

Pour chaque équation k :

```
Rabs_k = sum_C |R_C,k|                         [kg/s, N, N, W]
Rnorm_k = Rabs_k*tref / (sum_C V_C * Qref_k)
Rmax_k = max_C (|R_C,k|*tref / (V_C*Qref_k))
```

Les mêmes sommes L1 et maxima sont fournis pour chaque région, avec **son propre
volume**. Les critères imposent également `Rmax` : un grand Farfield ne peut donc
pas diluer un résidu d'ENGINE. Les résidus proviennent d'un nouvel assemblage de
l'état accepté, jamais seulement d'une variation de champ ou d'un résidu antérieur.
Les échelles sont physiques, pas les valeurs instantanées presque nulles d'une
composante. Elles doivent être adaptées au problème : une référence de gaz chaud
peut tolérer des variations relatives importantes dans un quasi-vide de masse
négligeable. Les diagnostics Mach/sections restent indispensables.

Le bilan stationnaire signé est `boundary_outflow_k-source_k`, source nulle sauf
pour le momentum radial. Il est normalisé par `Qref_k*Uref*Lref^2` : unités de
débit, force, force et puissance. Les variations des champs sont contrôlées en L1
volumique et en maximum par cellule, divisées par `Qref`. Les variations des deux
débits de bord utilisent `rho_ref*Uref*Lref^2`. L'imbalance relative de masse est
également fournie ; elle vaut zéro lorsque les deux débits sont exactement nuls.
Ce dernier rapport n'est pas utilisé seul pour un problème à débit quasi nul.

Une convergence exige simultanément, sur toutes les équations :

- résidus normalisés globaux, régionaux et maxima cellulaires sous la tolérance ;
- réduction suffisante des maxima initiaux s'ils étaient au-dessus de la tolérance ;
- variations globales et maximales des champs suffisamment faibles ;
- stabilité des débits et quatre bilans flux/source suffisamment faibles ;
- minimum d'itérations et plusieurs observations distinctes consécutives satisfaisantes.

Pour une équation déjà initialement sous la tolérance, aucune réduction impossible
d'un zéro ou d'un bruit d'arrondi n'est exigée : elle doit rester sous la tolérance.
Le diagnostic de réduction est `Rmax/max(Rmax_initial,tolerance)` ; ce dénominateur
explicite concerne uniquement ce diagnostic, pas le résidu ni le CFL.

`observe()` est idempotent à une même itération : des appels répétés ne fabriquent
pas d'observations stables. `observation_iteration` identifie l'état des diagnostics
entre deux contrôles planifiés. Une observation terminale ne reprend pas le solveur.

Statuts : `CONVERGED`, `MAX_ITERATIONS`, `WALL_LIMIT`, `STAGNATED`, `DIVERGED`,
`NUMERICAL_FAILURE`, outre `RUNNING`. Ce type distinct ne modifie pas les valeurs
sérialisées des anciens statuts du moniteur transitoire. `failure_message()`
conserve la cause d'un échec ; les états acceptés restent disponibles à l'export.

La divergence compare le pire résidu normalisé au pire résidu initial (ou à la
tolérance si initialement faible), avec un facteur de croissance explicite. La
stagnation utilise les **pics de fenêtres complètes, séparément pour les quatre
équations**. La première fenêtre fixe leurs références. Une baisse strictement
supérieure à `minimum_residual_progress` sur une équation encore non résolue
renouvelle sa référence et remet le compteur commun à zéro. Sinon, la référence
reste fixe afin que de petites baisses puissent s'accumuler. Une équation déjà
sous la tolérance ne peut pas masquer un plateau des autres ; si elle redevient
non résolue, sa référence est réinitialisée sans compter cette hausse comme un
progrès. L'arrêt exige `stagnation_windows` fenêtres consécutives sans progrès
significatif sur aucune équation non résolue, avec un résidu encore supérieur à
la tolérance. Par défaut, trois fenêtres sont requises : un ou deux pics
croissants peuvent appartenir à une oscillation acoustique amortie, comme sur
la tuyère 32x4. `stagnation_checks` fixe le nombre
d'observations par fenêtre, pas une limite d'itérations. Cette garde heuristique
reste configurable et peut arrêter une progression lente ; `STAGNATED` ne prouve
pas l'absence d'une racine. Les critères `CONVERGED` restent indépendants et inchangés.
La limite murale ne peut interrompre une opération élémentaire déjà engagée.

Défauts V1 : CFL `.05`, maximum `10000` itérations, `20` tentatives par itération,
`100` rejets cumulés, observation tous les `20`, minimum `20`, `3` contrôles
consécutifs ; tolérances résidu/champ/débit/bilan `1e-7`, réduction `1e-2`, limite
murale `60 s`. La stagnation utilise `50` observations par fenêtre, `3` fenêtres
sans progrès et une baisse significative de `1 %`. Le CFL par défaut reste `.05`.
Les gardes de divergence/stagnation sont documentées dans les options.
Aucune de ces valeurs n'est une garantie universelle pour RPA.

## Conservation : ne pas réutiliser le ledger physique

Chaque flux interne est toujours soustrait à l'owner et ajouté au neighbor.
Mais avec des pas locaux différents :

```
sum_C V_C*(Unew_C-Uold_C) = sum_C alpha*dtau_C*R_C
```

Cette somme n'est **pas** un bilan transitoire à un temps physique commun.
La modification des intégrales pendant les itérations n'est pas une fuite
physique mesurable par l'ancien ledger. Le contrôle stationnaire porte sur les
flux/source et le résidu non pondéré par les pas. Les tests vérifient précisément
cette distinction sur deux cellules de tailles différentes.

## Utilisation et exports

Avec `BUILD_VERIFICATION=ON` :

```sh
cmake --build build-release --target cfd_axisymmetric_euler_steady --parallel 2
./build-release/cfd_axisymmetric_euler_steady steady tube /tmp/steady-tube .05 12000 20 5
./build-release/cfd_axisymmetric_euler_steady transient tube /tmp/global-tube .05 12000 20 5
./build-release/cfd_axisymmetric_euler_steady steady nozzle /tmp/steady-nozzle .05 80000 20 20
./build-release/cfd_axisymmetric_euler_steady steady rpa /tmp/steady-rpa .05 50 5 5 /path/to/rpa-case
```

Arguments : `MODE CASE NEW_OUTPUT CFL MAX_ITER CHECK_INTERVAL WALL_SECONDS [RPA_CASE [uniform|quasi1d]]`.
`MODE=steady|transient`, `CASE=tube|nozzle|rpa`. Le pilote reprend les conventions
positionnelles des drivers Euler existants, pas un nouveau système de configuration.
Les critères avancés sont configurables via les options C++ ; aucun réglage du
chemin SIMPLE CLI n'est modifié. Le pilote borne ses invocations à 100000 itérations
et 120 secondes, mais **un budget CFD global externe reste nécessaire** pour une campagne.

`tube` est un mélange TRI/QUAD gradé avec perturbation de densité, entrée/sortie
supersoniques compatibles avec une solution uniforme. `nozzle` réutilise la petite
tuyère subsonique et sa projection quasi-1D existantes (16x2). Une quasi-1D n'est
pas une solution exacte 2D ; cette grille très grossière ne valide pas la précision
spatiale. `rpa` lit tous les paramètres du `meshDict` donné, sans modifier le cas.
Ses données nominales et hypothèses ambiantes exploratoires sont celles des essais
antérieurs (gaz équivalent, p_ext=1e-4 Pa, T_ext=300 K, vitesse extérieure nulle).

Pour RPA, le dernier argument facultatif sélectionne `uniform` (défaut historique,
identique bit à bit) ou `quasi1d` (choix explicite). Ce choix est enregistré dans
`summary.txt`. Exemple, à ne pas confondre
avec une certification de convergence :

```sh
./build-release/cfd_axisymmetric_euler_steady steady rpa /tmp/f4-quasi1d .05 100000 20 120 /chemin/vers/cas-F4 quasi1d
```

La préparation quasi-1D utilise le profil polygonal `inner_wall` du Mesh, son
rayon minimal et les conditions totales prescrites : relation aire–Mach,
branche subsonique avant le col et supersonique après, enthalpie totale et
isentropie constantes. Les valeurs ponctuelles aux centres cellulaires ont
une vitesse radiale nulle ; elles ne sont ni des moyennes cellulaires exactes,
ni une solution stationnaire 2D. Le débit critique est déduit de P0/T0/A*,
sans prescription indépendante. Seul ENGINE est remplacé ; les quatre autres
régions conservent exactement leurs conservatives historiques. Aucun choix
n'est transmis au solveur, aucun BC/CFL/critère n'est modifié. La qualité de
cette initialisation ne démontre pas la convergence asymptotique de F4.

VTU initial, intermédiaires espacés, final ; `snapshots.csv` indexe des **itérations**,
sans série PVD de faux temps physiques. `iterations.csv` fournit résidus, maxima,
bilans, extrema, pas, rejets et sections ; `regions.csv` identifie notamment ENGINE.
Le débit de section utilise le flux numérique existant. Dans le tube, `mid` est
une section de diagnostic, pas un col. Mach/p/T sont des moyennes des cellules
côté ENGINE, pas les traces exactes du Riemann de bord. Aucun fichier existant
n'est écrasé. Le CPU cumulé aux observations et le temps mural du stepping sont
enregistrés séparément ; leur différence n'est pas un coût exclusif d'une fonction.

Le mode `transient` de ce pilote est une référence physique bornée. Aux observations,
un inspecteur stationnaire temporaire réutilise l'assembleur pour obtenir exactement
les mêmes normalisations et bilans ; son coût est séparé du stepping. Il ne change
pas l'état transitoire. `MAX_ITERATIONS` n'est pas une certification stationnaire.
Un retour CLI zéro signifie une exécution contrôlée (limites et stagnation comprises),
pas nécessairement `CONVERGED` : consulter `summary.txt`. Les échecs numériques et
divergences renvoient 2, les erreurs de configuration ou d'I/O renvoient 1.
Comparer la réduction des mêmes résidus, pas des champs à un pseudo-temps égal à
un temps physique. Les itérations locales et physiques ne sont pas interchangeables.

Un champ transitoire accepté peut initialiser le solveur local par copie explicite.
Ce n'est pas une reprise complète du moniteur stationnaire. Les checkpoints
transitoires, leur format et leur lecture/écriture sont inchangés.

## Limites et validation

Les tests couvrent TRI/QUAD, axe/source radiale, vrai résidu, Jacobi, conservation
des flux internes, rejet après écriture partielle du candidat, seuils/options,
arrêts, protection contre la dilution Farfield et convergence d'un tube perturbé.
Un test protège les oscillations amorties contre une fausse détection de plateau.
La trajectoire transitoire est aussi comparée avant/après sur huit pas TRI/QUAD/mixte.

Ni positivité ni quelques observations décroissantes ne prouvent la convergence
sur F4. Limites : Euler explicite premier ordre, dépendance au CFL, réflexions des
frontières ouvertes, oscillations, stagnation et choix des échelles. Une détente
spatiale vers le quasi-vide reste un modèle continu exploratoire, sans validation
thermochimique, visqueuse ou raréfiée. Ne pas appeler « débit moteur stationnaire »
la vidange d'un état initial encore en évolution.
