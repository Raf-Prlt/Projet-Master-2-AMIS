<h1>PROJET M2 AMIS GROUPE A</h1>
Courtillat-Piazza Ludmila<br>
Guigal Alexis<br>
Hounguevou Thomas<br> 
Nzamboue Tisseu Hilary<br>
Préault Rafael<br>

Pour installer la libraire nauty: ouvrir un terminal dans le dossier NAUTY_2_8_8 et lancer la commande `make`.

Pour utiliser le programme, il faut lancer la commande `make` dans le dossier Transformation_Graphe.<br>
<br>
ATTENTION : le programme s'éxecute en fonction des informations présentent dans le fichier 'parametres.txt'.<br>
<br>
La première ligne correspond au mode d'éxecution : 

- mode 1 : le programme calcul les 35678 graphes de cycles, puis pour chacun, les 100 meilleurs scores de similarité et les écrit dans des fichiers, puis calcul les classes d'équivalences.

- mode 2 : affiche le fichier des 100 meilleurs scores de similarité d'une molécule x

- mode 3 : donne le score de similarité entre une molécule x et une molécule y

la ligne 2 correspond à l'ID de la molécule x.<br>
la ligne 3 correspond à l'ID de la molécule y.