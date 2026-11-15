# Graded-window barrier certificates — logiciel v1.1.0

[English](README.md)

Cette version étend le logiciel de reproduction aux **14 suites F4/F5** : contrôles à densité fixée, contrôles à haute densité et audits entiers/numériques. Un nouveau rejeu sur Windows a donné **14/14 PASS**, sans suite sautée. Les contrôles auxiliaires F3 ont donné **37/37 PASS**. Ces nombres désignent des exécutions de programmes ; ils ne remplacent pas une relecture mathématique humaine.

Ce dépôt contient les certificats logiciels, les données de test, les journaux et les instructions de reproduction. Le PDF du manuscrit et sa source LaTeX ne sont pas distribués ici. Les tags et archives v1.0.0 et v1.1.0 correspondent à des distributions contenant uniquement les certificats et leur documentation.

## Reproduire

```powershell
git clone --branch v1.1.0 https://github.com/MFRADJ/graded-window-barrier-certificates.git
cd graded-window-barrier-certificates
py -3.14 -m venv .venv
$python = '.\.venv\Scripts\python.exe'
& $python -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Installation des dépendances échouée' }
& $python -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Dépendances incompatibles' }
foreach ($mode in @('verify', 'all', 'f3')) {
    & $python reproduce.py --mode $mode --timeout 3600
    if ($LASTEXITCODE -ne 0) { throw "Échec du mode $mode" }
}
```

L'environnement exécuté est CPython 3.14.3, Windows 64 bits, python-flint 0.9.0, FLINT 3.6.0 et mpmath 1.3.0. Les commandes Linux/macOS sont dans le [README anglais](README.md) ; ces plateformes ne sont pas revendiquées comme exécutées pour cette release.

`verify` contrôle les empreintes des sources et de tous les fichiers distribués. `all` exécute les 6 suites standard, les 4 suites Arb et les 4 audits lourds. `f3` exécute les 37 vérifications auxiliaires, distinctes des 14 suites. `standard` et `arb` permettent des rejeux partiels. La limite de temps est définie par programme. Utiliser Python sans `-O` ni `PYTHONOPTIMIZE`.

Chaque exécution crée un nouveau sous-dossier de `results/local_runs/` avec ses sources copiées, journaux, artefacts, `CERTIFICATION_RESULT.json` et `SHA256_EXECUTION.txt`. Les sources distribuées et les journaux historiques ne sont pas réécrits. Une dépendance absente produit `SKIPPED` et un code non nul. Une limite de temps dépassée produit `TIMEOUT` et un code non nul. Une empreinte incorrecte interrompt le rejeu.

## Fichiers et portée

- Les quatre programmes historiques restent à leurs chemins d'origine : `certify_arb_independent.py` et les trois fichiers de `supplement/`.
- `certificates/high_density/` et `certificates/review/` contiennent les programmes supplémentaires, avec leurs octets d'origine.
- [La matrice des 14 suites](certificates/CERTIFICATE_MATRIX_14.csv) décrit le rôle et la méthode de chaque programme.
- [Les résultats v1.1.0](results/v1.1.0_2026-11-15/README.md) documentent les nouveaux rejeux.
- `results/f5_imported_2026-11-15/` conserve le résultat F5 initial : **6 PASS, 8 SKIPPED**, sans le présenter comme un succès complet.
- [La provenance](docs/PROVENANCE_v1.1.0.md) décrit la sélection depuis les archives, les licences et la protection de la v1.0.0.

Certaines suites associent coefficients entiers exacts et comparaisons logarithmiques flottantes. S2-bis comporte 31 contrôles Arb et 3 contrôles numériques. Les suites ne constituent donc pas toutes une preuve rigoureuse par intervalles. Le script heuristique `fixed_density_limit_V2-S4.py` est hors des 14 suites.

Le logiciel ne prouve ni ne réfute ETH, ne certifie pas formellement l'ensemble des théorèmes et ne fournit pas les opérateurs algébriques manquants. Le texte mathématique reste séparé du dépôt. Le logiciel et sa documentation sont sous [licence MIT](LICENSE). Utiliser [CITATION.cff](CITATION.cff) pour citer la version logicielle.
