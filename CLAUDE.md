# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Projet

Dahu est un serveur JSON-RPC d'analyse de données exposé via Tango, développé à l'ESRF
(initialement pour ID02). Le cœur (`src/dahu/`) est un moteur générique à plugins ; les
traitements scientifiques réels vivent dans `plugins/`, un sous-répertoire par ligne de
lumière (bm29, id02, id31, id15, id27...).

## Build et exécution

Le build utilise **meson-python** (pas setuptools). Le paquet n'est **pas importable depuis
les sources** : `src/dahu/__init__.py` lève une `RuntimeError` si `dahu.version` est absent
(ce module est généré à l'installation depuis `version.py`).

```bash
# Construire + lancer un script/commande contre le build local (patch de PYTHONPATH)
./bootstrap.py dahu-reprocess job.inp
./bootstrap.py ipython

# Installation
pip install .                # ou: pip wheel . && pip install dahu-*.whl

# Documentation Sphinx -> build/sphinx
./build-doc.py

# Paquet Debian
./build-deb.sh --help
```

`bootstrap.build_project()` fait `meson setup build` + `meson configure --prefix /` +
`meson install --destdir .`, puis devine le répertoire installé sous `build/lib/` : il vise
`build/lib/python3.X/site-packages` mais, si ce chemin n'existe pas, il remonte et redescend
par `os.listdir(...)[0]` — sur Debian on obtient en fait `build/lib/python3/dist-packages`.
C'est ce répertoire qui est ajouté au PYTHONPATH (surchargeable via `BUILDPYTHONPATH`).
Si `build/` a été produit par une version de meson plus ancienne, le build échoue avec
« Build data file ... references functions or classes that don't exist » : relancer
`meson setup --reconfigure build` (ou supprimer `build/`).

**Ajouter un fichier Python ne suffit pas** : il faut le déclarer dans le `meson.build` du
répertoire correspondant (`src/dahu/meson.build`, `plugins/meson.build`,
`plugins/bm29/meson.build`, etc.) sinon il ne sera jamais installé ni visible du moteur.

## Tests

```bash
python run_tests.py                        # construit puis lance dahu.test.suite
python run_tests.py --installed            # teste la version installée
python run_tests.py dahu.test.test_job.suite   # un module de test
python run_tests.py dahu.test.test_job.TestJob.test_plugin   # un seul test
python run_tests.py -v -c                  # verbeux + couverture
```

Les tests sont en `unittest` avec des fonctions `suite()` explicites : un nouveau test doit
être ajouté à la `suite()` de son module, son module à `src/dahu/test/test_all.py` **et** à
`src/dahu/test/meson.build` (sinon il n'est pas installé et `--installed` ne le voit pas).
Seul le noyau est testé (`test_job`, `test_plugin`, `test_cache`, `test_factory`,
`test_server`) ; les plugins ne le sont pas. `test_server` importe `dahu.server` **sans
installation Tango** : `stub_pytango()` place un faux module `PyTango` dans `sys.modules`,
ce qui permet d'exercer la logique pure (restitution des plugins indisponibles, `serialize`).

Piège actuel : `test_factory.py` et `test_server.py` sont déjà référencés par `test_all.py`
et par `meson.build`, mais ne sont pas suivis par git — sur un clone frais, `test_all` casse
à l'import.

Le lint en CI est `flake8` (`.github/workflows/python-package.yml`) ; `ruff` a été utilisé
récemment pour des passes de nettoyage mais n'est pas configuré dans `pyproject.toml`.

## Architecture

Chaîne d'exécution : **Tango (DahuDS) → Job (thread) → Factory → Plugin**.

- `src/dahu/plugin.py` — `Plugin` : constructeur vide, puis `setup()` (validation des
  entrées), `process()` (travail), `teardown()` (sortie + logs, toujours appelé même en cas
  d'échec), `abort()`. Entrée et sortie sont de simples dicts JSON-sérialisables
  (`self.input` / `self.output`). Les noms de ces méthodes sont surchargeables via
  `DEFAULT_SET_UP`/`DEFAULT_PROCESS`/`DEFAULT_TEAR_DOWN`/`DEFAULT_ABORT`.
  `plugin_from_function(f)` fabrique et enregistre une classe de plugin à partir d'une
  fonction sans état.
- `src/dahu/factory.py` — `Factory` (instance globale `plugin_factory`) : registre
  `{fqn: classe}` alimenté par le décorateur `@register`. Les noms de plugins sont
  **toujours en minuscules et pleinement qualifiés** (`module.classe`, ex. `example.cube`,
  `bm29.integratemultiframe`). Le chargement se fait par fichier (`importlib` sans polluer
  `sys.modules`), en cherchant dans : `$DAHU_PLUGINS` (séparé par `os.pathsep`), puis les
  répertoires passés au constructeur, puis `<site-packages>/dahu/plugins` — d'où
  l'installation de `plugins/` sous `dahu/plugins`. En pratique seuls `$DAHU_PLUGINS` et le
  répertoire par défaut comptent : `plugin_factory` est un singleton créé à l'import sans
  `plugin_path`, et la propriété Tango `plugins_directory` déclarée dans `DahuDSClass` n'est
  reliée à rien. Le gestionnaire de contexte `optional_plugin(fqn)` isole l'import et
  l'enregistrement d'un plugin : en cas d'échec, lui seul est désactivé et la raison est
  consignée dans `Factory.unavailable`, que `listPlugins` et `initPlugin` restituent.
- `src/dahu/job.py` — `Job(Thread)` : un job = un thread = une exécution de plugin.
  Identifiants entiers croissants, registre de classe `_dictJobs`, états
  `uninitialized/starting/running/success/failure/aborted`. `start()` instancie le plugin
  *avant* de démarrer le thread (si la fabrique renvoie `None` ou lève, le job passe
  directement en `failure` et les callbacks sont déclenchés — aucun thread n'est lancé).
  `abort_job_from_id(jobId)` (commande Tango `abort`) ne fonctionne que sur un job
  `running` : il bascule l'état en `aborted` puis appelle `plugin.abort()`, qui ne fait que
  positionner `self.is_aborted` — **c'est au plugin de tester ce drapeau** dans ses boucles,
  rien n'interrompt le thread. `clean()` sérialise entrée et sortie sur disque
  (`<workdir>/<jobid//1000>/<jobid>_<plugin>.inp|.out`, via `NumpyEncoder`) puis libère le
  plugin de la mémoire — après quoi les données ne sont plus lisibles que depuis le disque
  (`data_on_disk`). Les accesseurs statiques (`getDataOutputFromId`, `getStatusFromID`...)
  gèrent les deux cas. L'API historique est en camelCase avec une casse instable
  (`getJobFromID`/`getJobFromId`) : en renommant vers la convention PEP8, garder l'ancien
  nom — et ses variantes de casse — en alias de classe, comme pour `clean_job_from_id`.
- `src/dahu/server.py` — `DahuDS` / `DahuDSClass` : device Tango. Deux threads dédiés
  (`process_job` consomme `job_queue`, `process_event` pousse les events Tango
  `jobSuccess`/`jobFailure`). Commandes : `startJob([plugin, json])`, `waitJob`,
  `getJobState`, `getJobOutput`, `getJobInput`, `getJobError`, `abort(jobId)`, `listPlugins`,
  `initPlugin`, `cleanJob`, `collectStatistics`, `getStatistics`. L'attribut `serialize`
  force l'exécution séquentielle. `listPlugins` et `initPlugin` affichent aussi les plugins
  désactivés et la raison, lue dans `Factory.unavailable`.
- `src/dahu/cache.py` — `DataCache` : dict à taille bornée, **Borg par défaut** (toutes les
  instances partagent l'état) ; utilisé pour partager des objets coûteux entre plugins
  (ex. les intégrateurs azimutaux pyFAI dans `plugins/bm29/common.py`).
- `src/dahu/utils.py` — `get_workdir()` (répertoire de travail global `dahu_<isotime>`, créé
  une fois et mémorisé dans une globale), `get_isotime()`, `NumpyEncoder`,
  `fully_qualified_name()`.

Synchronisation entre plugins : `Plugin.wait_for(job_id)` attend un autre job (attribut de
classe `TIMEOUT`, 10 s par défaut) et échoue si celui-ci ne termine pas en succès.

### Points d'entrée

- `dahu-server` (alias historique `dahu_server`) → `dahu.app.tango_server:main`
- `dahu-register` → `dahu.app.tango_register:main` (enregistre le device dans la base Tango)
- `dahu-reprocess` → `dahu.app.reprocess:main` (rejoue des fichiers `.inp` hors Tango ;
  respecte `Plugin.REPROCESS_IGNORE` pour retirer des clés d'entrée non rejouables)

```bash
TANGO_HOST=localhost:10000 dahu-register --instance dahu
TANGO_HOST=localhost:10000 dahu-server dahu
dahu-server dahu -ORBendPoint giop:tcp::10001 -nodb -dlist id00/dahu/1 -v4   # sans base
```

`scripts/dahu_start` / `dahu_stop` sont les wrappers de production (logs sous `~/log/`, PID
dans `~/.dahu/pid`).

## `plugins/` — code d'analyse spécifique aux lignes de lumière

Le noyau ne sait rien de la physique : tout est ici. Un sous-répertoire (ou un fichier
`.py` à plat) par ligne, sans dépendance croisée. Ces fichiers ne sont jamais importés
comme `dahu.plugins.xxx` : la `Factory` les charge par chemin de fichier, donc le nom de
module vu à l'exécution est `bm29`, `id02`, `id27`... (pas `dahu.plugins.bm29`) — d'où les
FQN courts `bm29.integratemultiframe`, `id02.singledetector`.

| Famille | Ligne / domaine | Plugins exposés |
| --- | --- | --- |
| `bm29/` | BioSAXS (SAXS de protéines en solution) | `bm29.integratemultiframe`, `bm29.subtractbuffer`, `bm29.hplc`, `bm29.mesh` |
| `id02/` | TruSAXS (SAXS/WAXS résolu en temps) | `id02.metadata`, `id02.singledetector`, `id02.xpcs` |
| `id27.py` | diffraction haute pression (cellule diamant) | `id27.crysalisconversion`, `...fscannd`, `id27.xdiconversion`, `id27.average`, `id27.xdsconversion`, `id27.diffmap`, `id27.xdsprocessing` |
| `id31/` | diffraction de poudre | `id31.integrate`, `id31.integrate_simple` |
| `id15.py`, `id15v2.py` | diffraction/PDF | `id15.integratemanyframes`, `id15v2.integratemanyframes` (même nom de classe, FQN distincts) |
| `example.py` | référence/tests du noyau | `example.cube`, `example.square`, `example.noop`, `example.sleep` |
| `pyfai.py`, `focus.py` | démos anciennes | voir les pièges plus bas |

### Chaînages de jobs

Les plugins ne s'appellent pas entre eux : le séquenceur de la ligne (BLISS) soumet des
jobs indépendants et passe les identifiants des précédents dans la clé `wait_for`, que le
plugin consomme via `Plugin.wait_for(job_id)` dans son `setup()`. C'est le seul mécanisme
de synchronisation, et il échoue si le job attendu ne termine pas en `success`.

- **BM29** : `integratemultiframe` (une par acquisition, intègre N frames, CorMap pour
  détecter les frames équivalentes, moyenne) → `subtractbuffer` (moyenne les tampons
  équivalents, soustrait, Guinier/BIFT via `freesas`) → `hplc` (reconstruit le
  chromatogramme, NMF via `sklearn`) ou `mesh` (cartographie 2D).
- **ID02** : `metadata` (lit les compteurs du multiplexeur C216 via Tango, écrit le HDF5 de
  métadonnées) → `singledetector` (attend ce job via `metadata_job`, applique la chaîne
  dark/flat/solid-angle/distortion/normalisation puis l'intégration, cf. la clé `to_save`
  qui liste les étapes à écrire : `raw sub flat solid dist norm azim ave`).
- **ID27** : conversions et traitements largement délégués à des exécutables externes
  (`pyFAI-average`, `hdf2neggia`, `xds_par`, scripts `rsync`) lancés par `subprocess.run`.

### Conventions transverses

- **Chemins** : la convention ESRF `RAW_DATA` → `PROCESSED_DATA` (ancienne variante :
  `.../processed/`) est appliquée par substitution de chaîne pour déduire le répertoire de
  sortie quand `output_file` n'est pas fourni. Chaque plugin crée un sous-répertoire
  `gallery/` à côté de sa sortie ; `bm29/icat.py` en **déduit** proposition, ligne,
  échantillon et dataset par découpage du chemin — une arborescence non conforme casse
  silencieusement l'archivage.
- **Sorties** : un fichier HDF5/NeXus par job, écrit via la classe `Nexus`
  (`bm29/nexus.py`, `id02/nexus.py` — copies locales de celle de pyFAI/silx). L'entrée JSON
  complète est recopiée dans un `NXnote` « configuration », et chaque étape devient un
  `NXprocess` numéroté (`0_measurement`, `1_integration`, ...) avec un `NXdata` par défaut
  et un attribut `SILX_style`. Le dict de sortie du job ne contient que des chemins et
  quelques scalaires.
- **Archivage** : `to_pyarch` (fichiers `.dat`/PNG déposés dans pyarch pour ISPyB,
  `bm29/ispyb.py`, SOAP via `suds`), `send_icat()` (catalogue de données,
  `pyicat_plus`), `to_memcached()` (partage de résultats vers l'interface de ligne, port
  11211 en local). Les trois sont optionnels et ne doivent pas faire échouer le traitement.
- **Attente de fichiers** : les données arrivent souvent après la soumission du job ;
  `IntegrateMultiframe.wait_file()` scrute apparition **et** taille non nulle avec un
  `timeout` (clé d'entrée du même nom).
- **Intégrateurs pyFAI** : coûteux à construire, donc mis en cache (`DataCache`) et
  partagés entre jobs — cf. `bm29/common.py:get_integrator()` avec la clé
  `KeyCache(npt, unit, poni, mask, energy)`.

### Pièges vérifiés

- **Isolation des imports par famille** : `bm29/__init__.py` et `id02/__init__.py` ont été
  convertis au motif « un bloc `with optional_plugin(fqn):` par plugin », donc une dépendance
  manquante (`sklearn` pour `bm29.hplc`, `dynamix` pour `id02.xpcs`) ne désactive plus que le
  plugin concerné, la raison étant consignée dans `Factory.unavailable`. **Tout nouveau
  plugin doit suivre ce motif** ; `id31/__init__.py`, lui, importe encore tout directement
  (protections `try/except ImportError` + `logger.error` seulement sur pyFAI/fabio), et
  `id15*.py`/`id27.py` sont des modules à plat sans cette protection.
- **`DataCache` est un Borg par défaut** : `DataCache(3)` dans `id15.py`, `DataCache(6)`
  dans `id02/single_detector.py`, `DataCache(10)` dans `bm29/common.py` et `DataCache()`
  dans `id31/` sont **le même objet**, de taille fixée par le premier module chargé. Les
  clés ne se collisionnent pas (types différents), mais l'éviction est globale : charger
  des données ID02 peut évincer les intégrateurs BM29. Passer `borg=False` pour un cache
  réellement privé.
- `bm29/hplc.py` et `bm29/mesh.py` déclarent `def setup(self)` sans le paramètre `kwargs`
  de la classe de base. Cela passe parce que `Job._run_()` appelle la méthode sans
  argument, mais `plugin.setup({...})` en direct lève un `TypeError`.
- `plugins/pyfai.py` : `PluginIntegrate` et `PluginDistortion` ont leur `@register`
  commenté — seul `pyfai.integrate_simple` existe réellement ; `plugin_factory("pyfai.pluginintegrate")`
  renvoie `None`.
- `plugins/focus.py` définit `class Plugin(Plugin)` (masque la classe de base) et utilise
  `scipy.misc.imread`, supprimé de SciPy depuis la 1.3 : le plugin s'instancie mais son
  `process()` échoue. Code mort, à ne pas prendre comme modèle.
- `id27.py` lit `os.environ["HOME"]` **au moment de l'import** (chemin du greffon Neggia) :
  le module est inchargeable dans un environnement sans `HOME`. Il suppose aussi que les
  exécutables externes sont dans le répertoire de `sys.executable` (`PREFIX`).

## Écriture d'un plugin

```python
from dahu.plugin import Plugin, plugin_from_function
from dahu.factory import register

@register
class MonPlugin(Plugin):
    """Docstring obligatoire : la 1re ligne non vide est affichée par listPlugins."""
    def process(self):
        self.output["result"] = self.input.get("x", 0)
```

Conventions du dépôt :

- Un `plugins/<bl>/__init__.py` réenregistre ses classes sous des noms courts :
  `register(IntegrateMultiframe, fqn="bm29.integratemultiframe")`.
- Les dépendances scientifiques (pyFAI, fabio, h5py, hdf5plugin, freesas, sklearn, silx...)
  ne sont **pas** dans les `dependencies` du projet : elles s'importent dans le plugin, avec
  les précautions décrites ci-dessus.
- Utiliser `self.log_error(txt)` (lève une `RuntimeError`, le job passe en `failure`) et
  `self.log_warning(txt)` plutôt que `logger` seul : ces messages remontent dans
  `output["logging"]`.
- Un traitement long doit tester `self.is_aborted` régulièrement (et sortir proprement) :
  c'est le seul effet de la commande Tango `abort`.
- `do_profiling: true` dans l'entrée déclenche un profil cProfile écrit dans le workdir.
- Les sorties volumineuses vont dans des fichiers HDF5/NeXus (cf. `plugins/bm29/nexus.py`,
  `plugins/id02/nexus.py`) ; seul le chemin transite par le dict de sortie.
- En-tête de fichier standard : `__authors__`, `__contact__`, `__license__`,
  `__copyright__`, `__date__` (JJ/MM/AAAA), `__status__`, et `__version__` pour les plugins.
- Tester un plugin en développement sans réinstaller :
  `DAHU_PLUGINS=/chemin/vers/plugins ./bootstrap.py ...`

## Divers

`version.py` porte la version (schéma calendaire `MAJOR=2026, MINOR=3`) et est exécuté par
meson au moment du build. `GUI/`, `Lima_plugins/`, `datamodel/` et `example/` sont des
reliquats non construits et non installés, tout comme `src/dahu/plugin.py.orig` (suivi par
git mais absent de `meson.build`).
