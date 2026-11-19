# F4 : PLUME multibloc, génération et audit géométrique

Cette variante de **prétraitement** ne change ni le domaine physique, ni ENGINE,
ni Euler. Elle s'active explicitement ; les fichiers anciens restent utilisables.
Les cinq IDs `RpaRegion` et les huit groupes de frontières sont conservés.
Les sous-blocs sont des surfaces internes et des IDs de diagnostic, pas des BC.

## État de référence effectivement disponible

Le working tree fourni n'est pas le générateur horizontal archivé de la campagne
du 20 octobre. Les archives de 23990 cellules ne sont pas accessibles dans cet
environnement. Avec le DXF suivi, downstreamFactor=160, outerRadiusFactor=70,
bandThicknessFactor=1, plumeOutletRadiusFactor=18 et baseCellSize=.012, le chemin
historique actuel régénère **11495 cellules**, dont **2196 ENGINE**, et non 2184.
Il serait incorrect de présenter ces deux générateurs comme une référence
identique ou d'annoncer une accélération par rapport aux anciens calculs.

Le comparatif isolant cette modification utilise le chemin **plume finie déjà
présent avant la mission**, même DXF/domaine/résolution : plumeLengthFactor=40,
transitionWidthFactor=2, farfieldCellSize=.6, farfieldInterfaceCellSize=.15,
maximumGrowthRate=1.15. ENGINE et BAND sont identiques dans tous ces comparatifs.
Les archives CFD, cas privé et anciens checkpoints ne sont jamais ouverts en
écriture ni repris sur ce maillage. Une renumérotation liée à la nouvelle
topologie interdit de supposer une compatibilité des checkpoints.

## Topologie

Schéma non à l'échelle, pour une plume disposant de trois blocs :

```text
r_max  +---------------- FARFIELD supérieur (6) ----------+--------+
       |                                                |        |
       +------- enveloppe supérieure TRANSITION --------+  FF    |
       |  TRI (8), au-dessus de BAND et des blocs PLUME   | aval   |
       |   BAND (2)      / proche (3) / gradé (4) / loin(5)|  (7)  |
paroi  | ============== /-----------/---------/----------| cap TRI|
axe    + ENGINE (1) ----+-----------+---------+----------+--------+
                      sortie     x_near     x_cap      x_plume    x_aval
```

ENGINE demeure **une surface et une grille**, avec le même `wall_grid`, les
mêmes stations CAD/col, la même distribution radiale et les mêmes QUAD. BAND,
l'enveloppe et les deux blocs FARFIELD restent ceux du chemin fini existant.
PLUME devient un à trois blocs selon la place disponible. Une courbe verticale
commune est créée une fois et réutilisée par les deux surfaces adjacentes.
TRANSITION reste une seule surface TRI ; son contour supérieur PLUME devient
une chaîne de courbes conformes, sans nœud pendant ni traitement dans le solveur.

## Distributions et paramètres

| Clé `mesh` | Unité / défaut | Effet |
|---|---|---|
| `plumeNearLengthFactor` | rayon intérieur de sortie ; absent | Active le chemin multibloc ; longueur de la zone proche |
| `plumeAxialSpacingFactor` | adimensionnel ; 1 | Multiplie les cibles axiales héritées de la sortie et de la hauteur terminale moyenne |
| `plumeRadialGrowthRate` | adimensionnel ; absent | Rapport géométrique des intervalles radiaux sur la section terminale, entre 1 et maximumGrowthRate |

Il faut une plume finie HYBRID (`regions.plumeLengthFactor` et
`regions.transitionWidthFactor`). Les facteurs doivent être finis/positifs ;
les options multibloc ne sont pas silencieusement ignorées sur l'ancien chemin.
Les paramètres géométriques existants gardent leurs unités en rayons de sortie,
pas en diamètres. Il n'y a aucune coordonnée en mètres codée pour le moteur F4.

Avec h_e le dernier pas axial ENGINE, N_r le nombre de rangées PLUME, H_out sa
hauteur terminale et a=plumeAxialSpacingFactor :

* longueur proche = min(plumeNearLengthFactor * R_exit, longueur PLUME) ;
* cible axiale proche = a h_e ; au moins deux cellules, intervalle redistribué ;
* plafond axial terminal = max(a h_e, a H_out/N_r) ;
* au-delà : progression géométrique avec maximumGrowthRate, puis pas plafonné ;
* la fin du grossissement détermine l'interface du bloc lointain, sans position
  absolue arbitraire ;
* un reste plus court que deux cellules terminales est absorbé dans le bloc
  proche : aucun bloc minuscule n'est créé pour conserver trois noms ;
* le reste complet est redistribué pour atteindre l'extrémité exactement,
  sans minuscule dernière cellule. Le facteur est celui de la suite avant cette
  redistribution : ce n'est pas une garantie sur tous les sauts entre cellules.

À gauche, tous les nœuds ENGINE + lèvre + BAND sont conservés exactement.
N_r ne peut pas être choisi indépendamment à travers une interface tensorielle
conforme : il reste ENGINE + 1 intervalle de lèvre + BAND. L'option radiale
contrôle la **distribution**, pas le nombre. Par défaut le rapport historique
est plafonné par maximumGrowthRate sur ce nouveau chemin uniquement.
À chaque abscisse, les ordonnées sont l'interpolation linéaire entre les deux
sections. Les lignes correspondantes restent droites ; aucune forme précise
de jet n'est imposée comme paroi ni comme solution physique.

Exemple à ajouter au `mesh` d'une variante, pas au cas privé :

```text
plumeNearLengthFactor 4;  // deux diamètres de sortie
plumeAxialSpacingFactor 1;
// plumeRadialGrowthRate 1.05; // optionnel, sans modifier ENGINE
```

## Validation et diagnostics

La génération conserve le contrôle angulaire de chaque cellule et la preuve
QUAD mapped pour chaque cellule logique. Les tests comparent exactement les
coordonnées ENGINE/BAND avant/après, les nœuds de sortie/lèvre et les interfaces
PLUME partagées. Ils construisent Mesh et AxisymmetricGeometry, vérifient les
huit BC et la connexité. Aucun seuil existant n'est abaissé.

`cfd_rpa_hybrid_mesh` ne fait **aucun pas Euler**. Il produit un VTU global avec
`rpa_region`, `rpa_block`, volumes axisymétriques et tailles. Dans ParaView,
Threshold sur `rpa_block` isole ENGINE=1, BAND=2, PLUME=3..5,
FARFIELD=6/7 et TRANSITION=8. Les blocs absents n'ont pas de cellules.

```sh
cmake --build build-release --target cfd_rpa_hybrid_mesh cfd_rpa_nozzle_tests
timeout 90s ./build-release/cfd_rpa_hybrid_mesh VARIANTE NOUVEAU_DOSSIER
python3 -B verification/rpa_mesh_topology_report.py \
  NOUVEAU_DOSSIER/mesh.vtu NOUVEL_AUDIT --subsets
timeout 180s ctest --test-dir build-release \
  -R '^(rpa_nozzle|rpa_nozzle_integration|rpa_dxf_import|rpa_mesh_topology_report|gmsh_mesher)$' \
  --output-on-failure
```

L'audit indépendant produit JSON par région, CSV par cellule et vues SVG. Il
mesure min/max/moyenne et quantiles P1/P50/P99 : tailles, aires, volumes, qualité
de forme, qualité angulaire, aspect L²/A (TRI : L²/(2A)), non-orthogonalité et
rapport des racines d'aires des voisins. Pour l'angle, les faces de frontière
emploient centre-cellule vers centre-face ; le MeshReport standard ne compte
que les faces internes. Les formules sont donc explicitement distinguées.

`h_axi = V_axi / somme_f |S_axi,f|` localise les petites échelles géométriques.
**Ce n'est pas un pas CFL mesuré** : il manque les vitesses d'onde de l'état CFD.
Sur la référence courante, le minimum est près du col dans ENGINE et reste
inchangé (environ 1.7323e-4 m). Moins de cellules extérieures ne démontre donc
ni un plus grand pas temporel ni un gain de temps Euler.

## Domaine testé et limites

Tests synthétiques : profils conique, spline cubique, tuyère courte, longue et
rapport d'expansion différent ; h=.16/.08/.04 ; boîtes downstream=12/20/80,
rayon=8/10/30 ; plume courte, queue étroite, faible résolution, disparition des
blocs non réalisables ; rejets des options invalides et lecture meshDict.
Les dimensions chambre/profil sont portées par le **DXF** référencé dans
meshDict, pas par de nouveaux paramètres de chambre inventés.
F4 réel : h=.012/.008/.006, sans solve CFD. Cela ne couvre pas toutes les
géométries imaginables, notamment les profils non monotones en x ou invalides
déjà refusés par l'importeur, ni des rapports extrêmes non testés.

Les couches physiques optionnelles existantes sont conservées, sans ajouter
de couches limites Euler, de y+ imposé, d'injecteur ou de domaine solide. Les
futures physiques demanderont leur propre validation et éventuellement une
nouvelle résolution ; cette topologie ne les implémente pas. Le continuum dans
le quasi-vide n'est pas validé par le raffinement.

Résultats chiffrés, données avant/après, coûts et réserves de qualité : voir le
rapport de campagne de prétraitement livré avec les VTU. En particulier la
qualité de forme PLUME baisse localement lorsque l'espacement axial est plafonné
plus tôt ; ce compromis ne doit pas être dissimulé par la seule réduction de N.
