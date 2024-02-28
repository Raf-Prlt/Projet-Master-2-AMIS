<h1>PROJET M2 AMIS GROUPE A</h1>
Courtillat--Piazza Ludmila<br>
Guigal Alexis<br>
Hounguevou Thomas<br> 
Nzamboue Tisseu Hilary<br>
Préault Rafael<br>

Pour extraire les fichiers de la base de données: ouvrir un terminal dans le dossier TraitementDonnee et excécuter le fichier `main.py`

Pour installer la librairie nauty: ouvrir un terminal dans le dossier nauty2_8_8 et lancer la commande `make`.

Pour utiliser le programme, il faut lancer la commande `make` dans le dossier Transformation_Graphe.<br>
<br>
ATTENTION : le programme s'exécute en fonction des informations présentes dans le fichier 'parametres.txt'.<br>
<br>
La première ligne correspond au mode d'exécution: 

- **Mode 1** : Calcule les 35678 graphes de cycles, puis pour chacun, les 100 meilleurs scores de similarité, les écrit dans des fichiers 'similarite_id.txt' dans le dossier Résultats/similarite, calcule les classes d'équivalence et enregistre les résultats dans le fichier 'classesEquivalences.txt' dans le dossier Résultats.

- **Mode 2** : Affiche le contenu du fichier 'similarite_id.txt' contenant les 100 meilleurs scores de similarité d'une molécule $x$ et les molécules associées pour lesquelles ces scores ont été obtenus.

- **Mode 3** : Donne le score de similarité entre une molécule $x$ et une molécule $y$.

la ligne 2 correspond à l'ID de la molécule x.<br>
la ligne 3 correspond à l'ID de la molécule y.