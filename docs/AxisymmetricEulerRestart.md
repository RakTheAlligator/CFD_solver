# Euler F3d — checkpoint conservatif et étude CFL courte

**État au 18 octobre 2026.** Les résultats et checkpoints historiques sous
`/tmp` ont disparu. Les tableaux et commandes ci-dessous conservent le protocole
et les mesures rapportées ; ils n'impliquent pas que ces fichiers restent
disponibles. La nouvelle simulation repart explicitement de t=0, sans recréer
les comparaisons CFL perdues, et écrit hors dépôt sur un stockage persistant.

## Stockage persistant et nouveau calcul

Choisir un sous-dossier neuf de `$HOME/CFD_results/RPA/` sur un système de
fichiers persistant accessible en écriture ; ne jamais stocker ces données
dans Git ni écraser un ancien run. Conserver ensemble checkpoints, CSV,
VTU/PVD, métadonnées, statuts, rapports et provenance du binaire.

```sh
./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple "$HOME/CFD_results/RPA/f3d_20261018_cfl005_tau1" \
  300 0 0.05 200000 10000 0 0.000509251737436 250 20 600 \
  --physical-exports --checkpoint-every 10000
```

Cible autorisée : 1 t_a, sans prolongation automatique. Le flag
`--physical-exports` limite la visualisation à l'état initial et aux cinq
temps physiques planifiés. Une validation de reprise terminale doit utiliser
un autre dossier neuf, sans avancer au-delà de la cible. Les anciens repères
ne remplacent pas une comparaison de champs : leurs VTU ont disparu.

## Périmètre et revue

Voir [la revue F3](AxisymmetricEulerF3Review.md). Aucun changement des équations,
du Riemann Farfield, de Rusanov, des conditions physiques, du maillage, de SIMPLE
ou du cas privé. Aucun commit/push. Les fichiers existants et les résultats
F3c restent inchangés. Les modifications de sécurité concernent les contrats
des diagnostics, le plafond de rejets et la protection des sorties.

La revue préalable aux nouveaux objets de restart a passé **55/55 tests** en
Debug, Release et ASan+UBSan combinés. Le compte initial était 54 ; le test
Python du validateur puis le test C++ du checkpoint portent le total à 56.
Les validations finales passent **56/56 tests** en Debug, Release et
ASan+UBSan combinés ; clang-format, clang-tidy strict ciblé et `git diff --check`
passent également. Les logs sont conservés sous `/tmp/cfd-f3d-validated-*`.
`ASAN_OPTIONS=detect_leaks=0` : **LeakSanitizer n'est pas revendiqué**.

## Responsabilités et état sauvegardé

- `AxisymmetricEulerRestartState` possède les quatre tableaux conservatifs
  et les diagnostics du dernier pas accepté : temps, dt, CFL, compteurs,
  extrema, intégrales, flux de bord/source du pas et erreurs de bilan.
- Le constructeur de reprise du solveur valide cet état avec l'EOS et
  recalcule les diagnostics intégraux/extrema avant de l'accepter. Les buffers
  candidats, les résidus et sommes d'ondes sont reconstruits ; le choix du
  buffer courant n'a aucune influence mathématique. Les timers repartent de zéro.
- `EulerSteadyRestartState` conserve le champ de la précédente observation,
  **la vraie dernière observation**, le nombre de contrôles consécutifs et le
  marqueur d'observation. La reprise ne fabrique pas une observation pour
  remplir un historique et ne réinitialise pas les compteurs de convergence.
- `EulerCheckpointHistory` conserve les inventaires initiaux, les quatre
  intégrales temporelles des flux, la source radiale cumulée, les quatre
  échelles du ledger, le maximum d'erreur par pas et le planning d'exports.

Le cœur numérique ne connaît aucun chemin ni format de fichier. Les accès
`restart_state()` font des copies profondes **uniquement sur demande** ; aucun
buffer ou allocation de checkpoint par face/pas ordinaire. Le moniteur valide
ses diagnostics avant de modifier son champ historique, y compris en cas
d'exception. Sa mise à jour d'historique est une passe supplémentaire à chaque
observation, pas à chaque face ou pas temporel.

## Format version 1 et compatibilité

`cfd/io/EulerCheckpoint.hpp` fournit `write_euler_checkpoint()` et
`read_euler_checkpoint()`, séparés des constructeurs de restauration.

Texte versionné `CFD_EULER_CHECKPOINT 1`, locale classique, précision
`max_digits10<long double>` : au moins les 17 chiffres nécessaires aux doubles,
et 21 chiffres sur cette plateforme pour le ledger étendu. Les quatre
conservatifs sont écrits directement : **aucune reconstruction rho*u depuis VTU**.
Le lecteur vérifie aussi l'aller-retour canonique des valeurs `long double` ;
il refuse une précision étendue qui ne peut pas être reproduite sur sa plateforme.

La configuration stockée contient les coordonnées, types et ordres des cellules,
connectivités, faces owner/neighbor, groupes et noms de frontières, ainsi que
les métriques 2D effectivement utilisées. Elle contient gamma/R, CFL, politique
de retry/rejets, données de chaque BC, références/tolérances/cadence du moniteur
et contexte physique fourni par le driver (hypothèses ambiantes et régions RPA).
La comparaison de ces métadonnées est **complète et exacte**, sans `std::hash`.
Les budgets de temps/pas peuvent être étendus explicitement ; les paramètres
numériques et les critères de convergence ne peuvent pas changer silencieusement.

Une empreinte **FNV-1a 64 bits** couvre tous les octets du payload. Elle détecte
les corruptions accidentelles, **pas une modification malveillante** ; elle
ne sert pas de substitut à la comparaison des maillages/configurations.
Sont refusés notamment version inconnue, fichier tronqué, données supplémentaires,
checksum incorrect, cardinalité incorrecte, non-finis, gaz inadmissible et
historique de convergence incohérent. Les lectures sont bornées en taille
relativement au Mesh attendu avant allocation du payload.

L'écriture valide le snapshot, ferme/vérifie un fichier temporaire dans un
répertoire voisin exclusivement acquis, puis publie par `renameat2` avec
`RENAME_NOREPLACE` sur Linux. Aucun fichier existant n'est remplacé. Hors Linux,
la publication sans remplacement utilise un lien matériel ; une impossibilité
filesystem est une erreur explicite. Nettoyage temporaire au mieux lors d'une
exception. **Pas de garantie fsync contre une panne électrique**, pas de promesse
bit-à-bit entre compilateurs/arithmétiques flottantes différents.

Le reader conserve le ledger sans imposer un nouveau critère relatif universel
basé sur des flux nets. Un test radial donne un défaut axial absolu de
2.37e-21, mais un rapport net égal à 1 lorsque les forces s'annulent. Cette
normalisation n'est donc pas un critère universel d'admissibilité d'un checkpoint.
L'intégrité, les invariants de l'état et la continuation du ledger sont contrôlés
séparément ; les vérifications physiques restent celles du driver et des tests.

## Driver, sorties et arrêt contrôlé

Les arguments existants restent valides. Trois options nommées facultatives,
après les arguments positionnels, sont ajoutées :

```text
--checkpoint-every N     # cadence en compteur global de pas acceptés
--stop-after N           # arrêt contrôlé après N pas supplémentaires de ce segment
--restart FILE           # reprise dans un nouveau dossier externe vide
```

Lorsque le checkpointing est demandé, le dernier état est aussi sauvegardé
à la fin du segment. Noms : `checkpoint_STEP.chk`, ou `checkpoint_STEP_rREJECTS.chk`
si le compteur de rejets est non nul. Les anciens fichiers ne sont jamais
remplacés. Un échec d'export final ne doit pas empêcher la publication préalable
du checkpoint accepté. Le curseur de sortie désigne les futures cibles : un
segment repris exporte explicitement son état initial, même si une visualisation
avait échoué dans le segment précédent.

Chaque reprise possède ses propres CSV/VTU/PVD ; temps et numéros de pas sont
globaux. `assumptions.txt` indique `restart_from`, `segment_initial_time` et
`segment_initial_step`. Le PVD ne référence que les snapshots de son dossier.
Le snapshot initial de segment peut reproduire le dernier d'un segment précédent,
mais il est identifié comme tel, sans fusion implicite des chronologies.

Un état terminal Converged reste terminal. Une limite de temps/pas n'autorise
aucun pas supplémentaire sans budgets appropriés. Une extension du temps cible
définit de nouvelles cibles futures ; aucune sortie passée n'est réétiquetée.
La reproductibilité bit-à-bit exige les **mêmes contraintes temporelles futures**.

Aux frontières d'un segment, une observation non planifiée n'est pas injectée :
le VTU peut donc exister sans ligne de diagnostic à cet instant précis. Le
validateur reconnaît ce cas et contrôle la chronologie et le bilan du segment.
Le ledger du checkpoint demeure cumulatif depuis le vrai état initial.

## Tests de continuation et erreurs

- Petite onde synthétique : **100 pas = 50 + fichier + 50**, conservatifs
  identiques bit-à-bit, primitives, diagnostics temporels, résidus, planning,
  historique du moniteur et ledger étendu identiques. Le checkpoint au pas 50
  conserve l'observation réelle du pas 49, non une observation artificielle.
- Corruption, troncature, trailing content, version inconnue avec checksum
  recalculé, cardinalité, mauvais maillage de même taille, ordre des cellules,
  EOS, CFL, BC, hypothèses physiques, tolérance du moniteur, non-finis,
  état inadmissible, chemin inexistant, fichier existant, initial/terminal,
  vrais rejets, valeurs 1e-200 et zéro signé.
- RPA réel, 2240 QUAD : **200 = 100 + reprise + 100**, checkpoints et VTU finaux
  identiques octet par octet, cadence de diagnostic 17 (checkpoint entre deux
  observations). Validateur indépendant passé pour les segments.
- Échec VTU forcé au pas 100 : checkpoint valide conservé, reprise jusqu'à
  200 identique à la référence. Relecture d'un état à sa limite de pas :
  aucun nouveau pas et aucun historique réinitialisé.

Sur ces courts runs RPA : checkpoint environ **1.49 Mo**, écriture/validation
environ **32–39 ms**, lecture/validation environ **25 ms**, mesures unitaires
indicatives. Les copies et le staging texte restent linéaires et sont propres
aux événements I/O, pas à la boucle FV.

## Étude CFL : protocole et résultats

Même maillage, gaz, BC, initialisation et Rusanov ; seul CFL varie. Release
`-O3 -DNDEBUG`, un même binaire pour les quatre runs (empreinte conservée).
Cible **0.1 t_a = 50.9251737436 µs**, identique au premier palier F3c : évolution
non triviale déjà mesurable, sans recalcul de 20 minutes. Budget de boucle cumulé
mesuré **155.15 s**, inférieur à cinq minutes. Pas de balayage supplémentaire.
Un échantillon par CFL : tendances de coût, pas benchmark isolé ni optimum démontré.

| CFL | Pas | Rejets | CPU processus (s) | Boucle wall (s) | Mach col | Mach sortie | Q col (kg/s) | Q sortie (kg/s) | Engine (g) | Défaut normalisé max | Statut |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0.01 | 27629 | 0 | 88.94 | 85.512 | 0.470780 | 0.780705 | 0.290864271 | 30.7724454 | 8.78642105 | 3.17e-16 | TARGET_REACHED |
| 0.02 | 13816 | 0 | 47.16 | 43.725 | 0.470783 | 0.780712 | 0.290864123 | 30.7723589 | 8.78641763 | 5.05e-16 | TARGET_REACHED |
| 0.05 | 5528 | 0 | 21.01 | 17.368 | 0.470791 | 0.780736 | 0.290863677 | 30.7721015 | 8.78640779 | 1.12e-16 | TARGET_REACHED |
| 0.10 | 2766 | 0 | 12.32 | 8.548 | 0.470805 | 0.780774 | 0.290862919 | 30.7716768 | 8.78639216 | 1.15e-16 | TARGET_REACHED |

CPU inclut setup, import CAD et I/O ; wall ci-dessus désigne la boucle du driver.
Les coûts moyens sont 3.10 / 3.16 / 3.14 / 3.09 ms par pas. Les temps BC sont
inclus dans l'assemblage et les diagnostics Farfield dans le monitoring : ne
pas les additionner. Tous les états acceptés restent admissibles. Les cibles
sont atteintes, **aucun de ces calculs n'est stationnaire**.

| CFL | dt moyen (ns) | dt min naturel (ns) | dt max (ns) | p col (Pa) | p sortie (Pa) |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 1.843178 | 1.767002 | 1.946874 | 995651.073 | 501874.616 |
| 0.02 | 3.685956 | 3.534012 | 3.893787 | 995644.016 | 501871.104 |
| 0.05 | 9.212224 | 8.835103 | 9.734760 | 995622.832 | 501860.547 |
| 0.10 | 18.411126 | 17.670585 | 19.470501 | 995587.498 | 501842.909 |

Les pas tronqués aux exports ne sont pas utilisés pour le minimum naturel.
Masse globale finale 0.0106961919502827 kg et énergie 83918.0018368992 J,
identiques à l'arrondi près entre CFL ; cela ne démontre pas l'identité locale.
Les résidus changent peu : le maximum final reste autour de 8.915e-3.

Les conservatifs viennent des checkpoints, non des vitesses VTU. La comparaison
vérifie que seule la valeur CFL diffère dans les métadonnées numériques/physiques.
Pour chaque région G et champ q :

```text
L1_rel = sum_G(V*abs(q-qref)) / sum_G(V*abs(qref))
L2_rel = sqrt(sum_G(V*(q-qref)^2) / sum_G(V*qref^2))
```

Les normes absolues moyennes/RMS et Linf sont aussi fournies par le script.
Dénominateur nul : zéro si différence nulle, infini sinon ; aucune division
locale par une valeur proche de zéro. Les volumes sont reconstruits par
intégration polygonale, et Engine/extérieur sont normalisés séparément.

| Région | CFL | L2 rho | L2 rho ux | L2 rho ur | L2 rho E | L2 p | L2 Mach |
|---|---:|---:|---:|---:|---:|---:|---:|
| Engine | 0.02 | 3.85e-6 | 9.01e-6 | 3.80e-5 | 4.29e-6 | 4.51e-6 | 1.33e-5 |
| Engine | 0.05 | 1.54e-5 | 3.61e-5 | 1.52e-4 | 1.72e-5 | 1.80e-5 | 5.30e-5 |
| Engine | 0.10 | 3.47e-5 | 8.12e-5 | 3.42e-4 | 3.87e-5 | 4.06e-5 | 1.19e-4 |
| Extérieur | 0.02 | 1.46e-5 | 2.08e-5 | 2.48e-5 | 1.62e-5 | 1.47e-5 | 6.46e-4 |
| Extérieur | 0.05 | 5.82e-5 | 8.32e-5 | 9.89e-5 | 6.48e-5 | 5.85e-5 | 2.59e-3 |
| Extérieur | 0.10 | 1.31e-4 | 1.87e-4 | 2.22e-4 | 1.46e-4 | 1.31e-4 | 5.83e-3 |

Linf absolu Mach extérieur : 0.00526 / 0.02118 / 0.04804 pour CFL 0.02/0.05/0.10.
Les écarts croissent approximativement avec l'écart de pas, compatible avec
l'Euler temporel du premier ordre, sans constituer une preuve d'ordre absolu.
**CFL=0.01 n'est pas une référence temporellement convergée démontrée.**
La stabilité observée, la positivité et la conservation ne suffisent pas à
valider la précision ni la plume raréfiée.

## Prochaine expérience proposée — NON LANCÉE

Conserver CFL=0.01 pour les références historiques. Retenir **0.05 comme
candidat conservateur d'une nouvelle branche exploratoire**, pas comme optimum
universel : réduction d'environ cinq du nombre de pas, écarts régionaux
modérés sur cette fenêtre, marge par rapport au candidat 0.10.

Proposer la reprise du checkpoint 0.05 existant de 0.1 jusqu'à **0.5 t_a**,
pour comparer au palier 2 existant : environ 25000–35000 pas supplémentaires,
90–140 s de boucle estimées, plafond 240 s, 100000 pas globaux, 20 rejets globaux.
Diagnostics toutes les 250 étapes ; checkpoint tous les 10000 pas et en fin
de segment ; exports physiques restants 0.125/0.25/0.375/0.5 t_a, plus l'état
initial réel à 0.1 t_a. Budget estimé, non garantie de stabilité plus tardive.

Commande **à lancer seulement après accord**, vers un nouveau dossier :

```sh
./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/cfd-rpa-f3d-proposed-half-ta \
  300 0 0.05 100000 10000 0 0.000254625868718 250 20 240 \
  --physical-exports --checkpoint-every 10000 \
  --restart /tmp/cfd-f3d-cfl-005/checkpoint_5528.chk
```

Contrôler états admissibles, rejets, quatre bilans, normes régionales, débits,
masse Engine, distributions au col/sortie et régimes des trois Farfields.
Arrêt sur erreur numérique/I/O, budget atteint ou convergence réelle du
moniteur ; `TARGET_REACHED` reste distinct de Converged. Aucune poursuite
automatique vers 1 ou 2 t_a et aucune reconstruction exacte du palier 3 depuis VTU.
L'indépendance au domaine extérieur reste une étude séparée avant poussée/plume.

## Reproduction et conservation des artefacts

Étude courte (répertoire neuf à chaque CFL) :

```sh
./build-release/cfd_axisymmetric_euler_rpa_exploratory \
  cases/rpaNozzleSimple /tmp/NOUVEAU_DOSSIER \
  300 0 0.05 100000 10000 0 0.0000509251737436 250 20 150 \
  --physical-exports --checkpoint-every 100000
python3 verification/validate_rpa_euler_output.py /tmp/NOUVEAU_DOSSIER
python3 verification/compare_rpa_euler_cfl.py \
  /tmp/cfd-f3d-cfl-001 /tmp/cfd-f3d-cfl-002 \
  /tmp/cfd-f3d-cfl-005 /tmp/cfd-f3d-cfl-010
```

Résultats conservés sous `/tmp/cfd-f3d-*`, logs de CTest/build/tidy également.
`/tmp/cfd-f3d-cfl-comparison.txt` contient toutes les normes L1/L2/Linf.
Les anciennes données `/tmp/cfd-rpa-f3c-stage*` sont préservées par SHA256.
**/tmp peut être nettoyé par le système** : archiver les résultats et checkpoints
importants dans un emplacement persistant externe choisi par l'utilisateur.
Aucune copie ou suppression automatique des anciens résultats n'est faite.
