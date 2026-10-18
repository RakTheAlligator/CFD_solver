# Revue F3 et préparation F3d

Compte rendu historique de la revue. Les résultats temporaires qui avaient
alors été protégés par SHA256 ne sont plus disponibles au 18 octobre 2026.
Le stockage persistant et le protocole de relance sont décrits dans
[AxisymmetricEulerRestart.md](AxisymmetricEulerRestart.md).

Revue sur `feat/compressible-axisymmetric`, HEAD initial `f605229` ; le diff
F3a/b/c non committé préexistant est conservé. Aucun changement SIMPLE ni du
cas privé `cases/rpaNozzleSimple`. Les résultats F3c sont protégés par SHA256.

## Fichiers et interactions examinés

- `AxisymmetricGeometry.hpp/.cpp` : métriques, identités de surface/source.
- `ConservativeField.hpp/.cpp`, `CaloricallyPerfectGas.hpp/.cpp` : SoA,
  ownership, admissibilité, conversions et arithmétique élargie.
- `EulerFlux`, `EulerSlipWallFlux`, `EulerBoundaryConditions`,
  `EulerCharacteristicBoundaryFlux`, `EulerFarfieldBoundaryFlux` (headers et
  sources) : orientations, gaz unique, traces et vitesses d'onde.
- `AxisymmetricEulerSolver.hpp/.cpp` : assemblage, source, CFL, buffers,
  backtracking, commit et diagnostics.
- `EulerSteadyMonitor.hpp/.cpp`, `EulerSteadyCsvWriter.hpp/.cpp` : résidus
  courants, inventaires, historique et normalisations.
- `verification/support/EulerBoundaryWave.hpp`, `EulerRpaExploratory.hpp`,
  `EulerRpaTemporal.hpp`, les drivers Farfield/RPA et le validateur Python.
- Tests flux, paroi, caractéristiques, Farfield, solveur, Sod, expansion,
  moniteur, initialisation et sections ; CMake et configuration sanitizers.

F3a : Farfield, type BC, dispatch, tests et démonstration acoustique/backflow.
F3b : initialisation RPA, driver exploratoire, sorties et validateur.
F3c : sorties physiques, sections, masse Engine et étude temporelle.
Le répertoire privé est indépendant et n'est pas ajouté au suivi Git.

## Constats classés

### 1. Erreurs corrigées

- `RpaFaceSection` : une géométrie provenant d'un autre Mesh, ou un appel
  `sample()` avec un Mesh/champ/BC incompatible, pouvait indexer des données
  incorrectes. Identité Mesh et cardinalités désormais contrôlées ; tests
  avec un autre Mesh de même taille et un champ court.
- `AxisymmetricEulerSolver::step()` / driver RPA : le contrôle cumulatif
  après un pas accepté pouvait dépasser le budget de rejets pendant les
  tentatives. Option explicite `maximum_rejected_steps`, contrôlée avant
  toute nouvelle tentative ; le dernier état accepté reste intact.
  Le défaut générique illimité conserve les comportements antérieurs.
  Régression sur le cas radial qui exige effectivement un backtracking.
- `write_backflow_demo()` : les fichiers de démonstration pouvaient écraser
  un dossier existant. Un dossier non vide est maintenant refusé.
- Validateur RPA : contrôle explicite ajouté entre le temps final exporté
  et la dernière ligne de l'historique accepté. Le simple fait qu'un temps
  ait été accepté antérieurement ne suffit pas. Régression Python dédiée.

Aucune formule de flux, source, EOS ou mise à jour n'a changé. Les trajectoires
valides sans rejet restent inchangées ; aucun résultat antérieur n'est effacé.

### 2. Mauvaises pratiques / limites de responsabilité

- Le driver RPA est long et rassemble orchestration et diagnostics. Il est
  spécifique à cette expérience : aucune refactorisation générale entreprise.
- L'acquisition du dossier de résultats suppose un seul producteur ; elle
  n'est pas un protocole de coordination de processus concurrents.
- Un échec d'I/O peut laisser un CSV/VTU incomplet ; le statut d'échec et le
  validateur doivent être consultés. Un checkpoint valide ne doit pas dépendre
  de la publication réussie d'un fichier de visualisation.

### 3. Risques de performance (non optimisés)

- Le Riemann Farfield réévalue l'état extérieur et résout la pression étoile
  en logarithmes à chaque appel. Aucune approximation ni cache ajouté.
- Les diagnostics de frontière recalculent les traces à chaque observation.
  Leurs timers recouvrent le monitoring ; les BC de production sont déjà
  incluses dans l'assemblage. Ces temps ne doivent pas être additionnés.
- Flush de `steps.csv` à chaque pas : coût I/O potentiel, sans optimisation
  arbitraire dans ce lot. Aucun tableau alloué par face ou par pas valide.

### 4. Améliorations facultatives

- Étude distincte d'indépendance au domaine extérieur avant une validation
  quantitative du panache ou de la poussée.
- Archivage persistant externe des résultats importants : `/tmp` est
  temporaire. Aucun déplacement automatique des données effectué.

### 5. Préférences stylistiques

Pas de renommage massif, de framework ou de nouvelle hiérarchie pour harmoniser
les petits helpers de validation ; ce ne sont pas des erreurs scientifiques.

## Conclusions mathématiques et limites

Les signes correspondent à `V dU/dt = -sum(A Fhat) + [0,0,2*pi*A_cell*p,0]`.
Chaque flux interne est réutilisé avec des signes opposés. L'axe est ignoré
avant normalisation par une surface nulle. Le gaz évalue l'énergie interne
avant d'accepter un état ; les candidats rejetés ne sont pas publiés.

Rusanov est consistant et conservatif mais diffuse les contacts. Le flux
Farfield est le flux physique de la trace Riemann normale à xi=0 (pas une
seconde dissipation Rusanov) ; le côté du contact porte l'entropie et la
vitesse tangentielle. Les chocs et éventails soniques sont échantillonnés
selon leurs vitesses. Le vide généré au Farfield est explicitement refusé ;
la paroi traite séparément son vide physique. Ce n'est pas une BC exacte
multidimensionnelle non réfléchissante.

Le CFL de faces n'est pas une preuve de positivité avec source radiale :
le backtracking reste nécessaire. Les comparaisons de Sod et de l'expansion,
les références analytiques des flux et les intégrales polygonales du
validateur fournissent des contrôles différents ; ils ne constituent pas
une preuve formelle de tout le logiciel. Les limites de représentation et
de séparation énergie interne/cinétique à très grand Mach restent réelles.

Les Mach de section sont des estimations cellulaires, les débits sont des
flux numériques de faces. Résidu stationnaire, inventaire transitoire et
défaut du ledger sont distincts. Les hypothèses ambiantes T=300 K, vitesses
nulles et gaz équivalent restent exploratoires, non des données spatiales
mesurées. Aucun défaut mathématique bloquant identifié dans le périmètre revu.

Les validations avant restart et les validations finales F3d sont rapportées
dans `AxisymmetricEulerRestart.md` ; LeakSanitizer n'est pas revendiqué
lorsque `ASAN_OPTIONS=detect_leaks=0` est utilisé.
