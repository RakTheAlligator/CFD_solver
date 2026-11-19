# RPA F4 — plume finie, raccord TRI compact et Farfield mapped

F4 ne change ni Euler, ni l'EOS, ni les flux, ni les huit conditions aux limites.
Le cas privé et les résultats historiques restent en lecture seule. Utiliser une
**variante séparée** de `system/meshDict`, jamais des constantes dans le solveur.

## Géométrie physique, régions et discrétisation

La géométrie physique demeure le profil aplati importé, sa paroi extérieure,
la lèvre, l'axe et la boîte. Les facteurs historiques multiplient le **rayon
intérieur de sortie** Re ; De = 2 Re. `upstreamFactor auto` conserve l'amont x=0
sur le cas étudié : il ne signifie pas une extension devant le moteur.

La nouvelle topologie est activée par `regions.plumeLengthFactor` et
`regions.transitionWidthFactor`, à fournir ensemble. **Absents : chemin historique
inchangé.** Elle exige `cellType hybrid`.

ENGINE, BAND, PLUME sont des grilles mapped. Deux blocs mapped couvrent FARFIELD :
un bloc supérieur, et un bloc aval raccordé à son côté droit. TRANSITION est
l'unique surface libre Gmsh, exclusivement TRI. Elle passe au-dessus de BAND et
PLUME, puis derrière l'extrémité de PLUME jusqu'à l'axe. PLUME et FARFIELD n'ont
ainsi aucune interface commune, même en aval. Toutes les subdivisions sont
internes, sans nouvelle BC.

Notations : xe, Re désignent la sortie réelle, L la longueur de plume et δ la
largeur configurée, en mètres. Le bas du bloc supérieur vaut :

* avant xe : max(profil extérieur + épaisseur BAND) + δ ;
* entre xe et xe+L : interpolation linéaire vers rayon terminal PLUME + δ ;
* entre xe+L et xe+L+δ : rayon terminal PLUME + δ.

Le bloc aval commence à xe+L+δ. Ces interfaces de maillage ne sont pas des parois.
δ doit laisser une largeur aval et une hauteur supérieure strictement positives.
Ce découpage n'impose plus de contact entre le haut de PLUME et le bas de FARFIELD.
La borne Fmin de l'ancienne pointe commune ne s'applique pas à ce raccord avec
cap TRI : elle reste vérifiée dans la topologie historique. La validation angulaire
de **toutes les cellules générées** reste inchangée, sans abaissement de seuil.

## Paramètres

| Clé | Unité / défaut | Contrat |
|---|---|---|
| `regions.plumeLengthFactor` | Re ; absent | Longueur organisée après xe, strictement avant l'aval |
| `regions.transitionWidthFactor` | Re ; absent | Largeur δ positive au-dessus/derrière la plume ; obligatoire avec plume finie |
| `mesh.engineRadialCells` | nombre ; absent | Au moins 2, compatible avec les couches explicites ; sinon sélection par taille de base |
| `mesh.bandRadialCells` | nombre ; absent | Même validation pour BAND |
| `mesh.engineAxialSpacingFactor` | adimensionnel ; 1 | Multiplicateur positif du pas axial dérivé de la hauteur moyenne |
| `mesh.preserveWallVertices` | booléen ; false | Conserve chaque sommet du profil aplati comme station, sans modifier le CAD |
| `mesh.farfieldInterfaceCellSize` | m ; baseCellSize | Taille positive côté mapped de TRANSITION, au plus farfieldCellSize |
| `mesh.farfieldCellSize` | m ; baseCellSize | Plafond positif des espacements Farfield et du champ de taille TRI |
| `mesh.maximumGrowthRate` | adimensionnel ; 1.2 | Strictement >1 ; progression des lignes Farfield et pente TRI = facteur−1 |
| `mesh.plumeRefinementRadiusFactor` | Re ; absent | Enveloppe de raffinement TRI à xe, au plus le rayon terminal PLUME |
| `mesh.plumeRefinementCellSize` | m ; absent | Taille positive de cette enveloppe ; obligatoire avec son rayon |

Le pas axial du bloc supérieur est proche de farfieldInterfaceCellSize et inclut
les coins de l'enveloppe. Les hauteurs de rangées croissent géométriquement puis
sont plafonnées. Le bloc aval réutilise exactement la distribution du côté droit
du bloc supérieur ; ses colonnes croissent vers l'aval avec le même plafond.
Les suites sont redistribuées sur leur intervalle entier, sans minuscule dernier
élément. Dans le bloc supérieur, les hauteurs sont normalisées puis adaptées à
la hauteur de chaque colonne : le plafond est une cible sur la colonne de référence,
pas une borne stricte sur chaque arête physique. Ce choix simple accepte une
anisotropie mesurée dans les QUAD lointains.
maximumGrowthRate ne garantit pas le rapport des racines d'aires entre tous voisins,
notamment TRI/QUAD ou aux stations CAD obligatoires.

Les deux paramètres d'enveloppe TRI ne pilotent pas les grilles mapped. Les options
d'extérieur exigent une plume finie et les configurations incohérentes sont rejetées.
Avec rangées et tailles explicites, diminuer baseCellSize seul ne raffine pas ces
overrides. Les valeurs par défaut absentes préservent les cas historiques.

## Variante étudiée et choix mesurés

```text
FoamFile { format ascii; class dictionary; object meshDict; }
geometry
{
    type rpaNozzle;
    engine { file "constant/geometry/engine.dxf"; wallThicknessFactor 0.02; }
    box { upstreamFactor auto; downstreamFactor 160; outerRadiusFactor 70; }
}
regions
{
    bandThicknessFactor 0.05;
    plumeOutletRadiusFactor 18;
    plumeLengthFactor 40;
    transitionWidthFactor 0.9;
}
mesh
{
    cellType hybrid;
    baseCellSize 0.004;
    minCellQuality 0.05;
    engineRadialCells 36;
    bandRadialCells 8;
    engineAxialSpacingFactor 1.5;
    preserveWallVertices true;
    farfieldCellSize 0.30;
    farfieldInterfaceCellSize 0.05;
    plumeRefinementRadiusFactor 8;
    plumeRefinementCellSize 0.020;
    maximumGrowthRate 1.2;
}
```

Adapter le chemin DXF dans une variante extérieure au cas privé. Sur le profil
étudié, De=.285816514 m et xe=.650175702 m : aval=xe+80 De=23.515496822 m,
rayon extérieur=35 De=10.003577990 m, amont=0. La plume finit à xe+20 De,
avec rayon terminal 9 De (facteur historique 6 → 18, exactement ×3).
Son cœur QUAD débute sur les nœuds de sortie/lèvre/bande sans les déplacer.

La taille de base .004 m correspond à environ 36 rangées à la sortie. Le pas axial
1.5 fois la hauteur moyenne évite une sur-résolution uniforme, mais préserver tous
les sommets conduit encore à 514 colonnes ENGINE. C'est le coût de la fidélité au
profil aplati, pas un objectif de nombre de cellules. BAND est une zone de maillage,
pas une couche limite physique ni un modèle y+.

Recherche bornée de δ : la taille d'interface .05 m et la pente .2 donnent une
longueur de grossissement indicative .25 m. Les facteurs 1.8, .9, .45 et 3.6 ont été
comparés, soit des largeurs proches de cette longueur, de ses moitiés et de son double.
Les qualités de forme minimales TRI mesurées sont respectivement .559, .550, .480,
.570. Le facteur .9 est le plus compact testé atteignant l'objectif de conception
.5 ; **ce n'est pas un nouveau seuil de rejet utilisateur**. Transition occupe
1.216 m², environ .52 % de l'aire méridienne du domaine.

La taille Farfield .75 m de la première tentative donnait des QUAD mapped d'aspect
maximal ~21 ; .30 m abaisse le maximum global à ~9.9, pour 9510 QUAD Farfield.
Un essai BAND à 12 rangées réduit le saut maximal de tailles 4.60 → 3.72, mais
augmente le coût et réduit les minima de qualité PLUME/TRANSITION ; il n'est pas
retenu. Le candidat BAND8 garde une qualité TRI min/P1/médiane .550/.770/.984.
Le saut 4.60, localisé près du col, reste une limite explicitement mesurée.

## Preuves mapped, contrôles et exports

Les grilles sont des tableaux lignes/colonnes explicites, installés sur des surfaces
transfinite recombinées, pas une recombinaison libre de TRI. La génération vérifie
chaque QUAD contre les quatre coordonnées de sa cellule logique, puis exige TRI
uniquement dans TRANSITION. Les deux blocs Farfield partagent leurs nœuds et segments.
Le rapport donne les colonnes/rangées/cellules de chacun des cinq blocs.

Le préflight, sans solve, lit le meshDict réel :

```bash
cmake --build build-release --target cfd_rpa_hybrid_mesh
./build-release/cfd_rpa_hybrid_mesh VARIANT_CASE NEW_EMPTY_MESH_OUTPUT
```

Il produit mesh.vtu (régions, tailles, volumes axisymétriques) et mesh-report.txt.
Les contrôles indépendants des artefacts reconstruisent les blocs logiques à partir
des IDs réellement exportés, les faces partagées, les volumes et la connexité.
Les validateurs Python lisent les vrais types/offsets VTK, pas un pas fixe de 4.

Groupes conservés : 0 inlet (TotalConditionsInlet), 1 axis (Axis), 2 inner_wall,
3 outer_wall, 4 nozzle_lip (SlipWall), 5 farfield_upstream, 6 farfield_radial,
7 farfield_downstream (Farfield). Données physiques et hypothèses ambiantes :
voir AxisymmetricEulerRpaExploratory.md. Aucun débit imposé ni changement d'EOS.

Les tests utilisent des DXF synthétiques publics, jamais le cas privé. Les résultats,
vues et variantes F3 / premier F4 / nouveau F4 sont archivés séparément hors dépôt.
Le démarrage autorisé est limité à 500 pas, CFL=.05, sans ancien checkpoint.
Il ne démontre ni convergence stationnaire, ni indépendance spatiale/aux frontières,
ni validité du modèle continu pour le panache raréfié. Le coût long dépend du CFL,
de l'état de l'écoulement et des conditions CPU ; ne pas l'extrapoler sans réserve.
