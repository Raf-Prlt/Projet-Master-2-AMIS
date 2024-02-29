#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "comparaison.h"

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de cycles
double comparaisonNbreDeCycles(struct g_cycles g1, struct g_cycles g2){
    double degreDeSimilarite;
    int max;
    int min;
    if (g1.nb_cycles < g2.nb_cycles){
        min = g1.nb_cycles;
        max = g2.nb_cycles;
    }else{
        max = g1.nb_cycles;
        min = g2.nb_cycles;
    }
    degreDeSimilarite = (double) min / max;

    return degreDeSimilarite;
}

// Fonction pour comparer deux graphes de cycles sous le critère de la taille des cycles
 double comparaisonTailleDeCycles(struct g_cycles g1, struct g_cycles g2){
        double degreDeSimilarite;
        double cyclesIdentiques = 0.0;
        double max;
        if (g1.nb_cycles < g2.nb_cycles){
            max = g2.nb_cycles;
        }else{
            max = g1.nb_cycles;
        }
        int i = 0;
        int j = 0;
        
        while (j < g2.nb_cycles && i < g1.nb_cycles) {

          if (g1.generateur[i].taille == g2.generateur[j].taille){
            cyclesIdentiques++;
            j++;
            i++;
          } else if (g1.generateur[i].taille < g2.generateur[j].taille) {
            i++;
          } else if (g1.generateur[i].taille > g2.generateur[j].taille) {
            j++;
          }

        }
        degreDeSimilarite = (double) cyclesIdentiques / max;
        return degreDeSimilarite;  
}

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de voisins
double comparaisonNbreVoisinsDeCycle(struct g_cycles g1,struct g_cycles g2){
    double degreDeSimilarite;
    double max;
    double min;
    double degre1 = 0;
    double degre2 = 0;

    for (int i = 0; i < g1.nb_cycles; i++) {
        degre1 += g1.generateur[i].degre;
    }
    for (int i = 0; i < g2.nb_cycles; i++) {
        degre2 += g2.generateur[i].degre;
    }

    degre1 = (double) degre1 / g1.nb_cycles;
    degre2 = (double) degre2 / g2.nb_cycles;

    if (degre1 < degre2){
        min = degre1;
        max = degre2;
    }else{
        max = degre1;
        min = degre2;
    }

    if (max == 0) {
        degreDeSimilarite = 1;
    } else {
        degreDeSimilarite = (double) min / max;
    }
    return degreDeSimilarite;
}

// Fonction pour comparer deux graphes de cycles sous le critère du poids moyen des liaisons entre les cycles
double comparaisonPoidsDesArretes(struct g_cycles g1, struct g_cycles g2){
        double degreDeSimilarite;
        double poids1 = 0;
        double poids2 = 0;
        double max;
        double min;

        if (g1.nb_liaisons == 0 && g2.nb_liaisons == 0) {
            return 1;
        }

        if (g1.nb_liaisons == 0 || g2.nb_liaisons == 0) {
            return 0;
        }  

        for (int i = 0; i < g1.nb_liaisons; i++) {
            poids1 += g1.aretes[i].Poids;
        }

         for (int i = 0; i < g2.nb_liaisons; i++) {
            poids2 += g2.aretes[i].Poids;
        }

        poids1 = (double) poids1 / g1.nb_liaisons;
        poids2 = (double) poids2 / g2.nb_liaisons; 

        if (poids1 < poids2){
            min = poids1;
            max = poids2;
        }else{
            max = poids1;
            min = poids2;
        }
            
    if (max == 0) {
        degreDeSimilarite = 1;
    } else {
        degreDeSimilarite = (double) min / max;
    }
    return degreDeSimilarite;
}

double similarite(struct g_cycles g1, struct g_cycles g2){
    double degreDeSimilarite;
    double degre = 0;
    degre += comparaisonNbreDeCycles(g1, g2);
    degre += comparaisonTailleDeCycles(g1, g2);
    degre += comparaisonNbreVoisinsDeCycle(g1, g2);
    degre += comparaisonPoidsDesArretes(g1, g2);

    degreDeSimilarite = (double) degre / 4;
    return degreDeSimilarite;

}

// Calculer les similarités entre les graphes de cycles et stocker le résultat dans un fichier
void calculSimilarite(struct g_cycles **liste, int size){
    //Création du sous dossier
    #ifdef _WIN32
        system("mkdir Résultats 2> nul"); // Pour Windows
        system("mkdir Résultats/similarite 2> nul"); // Pour Windows
    #else
        system("mkdir -p Résultats 2> /dev/null"); // Pour Linux/Unix
        system("mkdir -p Résultats/similarite 2> /dev/null"); // Pour Linux/Unix
    #endif

    // Calcul du degré de similarité
    for (int i = 3; i < size; i++)
    {
        double *tab = malloc(100 * sizeof(double));
        //Définition d'une sous liste
        struct g_cycles **sliste = malloc(100 * sizeof(struct g_cycles *));

        for (int y = 0; y < 100; y++) {
            tab[y] = 0.0;
            sliste[y] = liste[y+3];
        }
        
        for (int j = 3; j < size; j++)
        {
            double val = similarite(*liste[i], *liste[j]);

            double min = 100.0;
            int indiceMin = -1;

            for (int t = 0; t < 100; t++) {

              if(min >= tab[t]) {
                min = tab[t];
                indiceMin = t;
              }

            }

            if (val >= min) {
              tab[indiceMin] = val;
              sliste[indiceMin] = liste[j];
            }
        }
        
        // Tri décroissant par ligne
        for (int k = 0; k < 100; k++)
        {
            for (int j = 0; j < 100; j++)
            { 
                if(tab[k] > tab[j])
                {
                    //struct g_cycles *x = malloc(sizeof(struct g_cycles));
                    struct g_cycles *x = sliste[k];
                    double c = tab[k];

                    tab[k] = tab[j];
                    sliste[k] = sliste[j];

                    tab[j] = c;
                    sliste[j] = x; 
                    //free(x);   
                }
            } 
        }

        for (int t = 0; t < 100; t++) {

          if (sliste[t]->Id == liste[i]->Id) {
            struct g_cycles *x = sliste[t];
            sliste[t] = sliste[0];
            sliste[0] = x; 
          }
        }

        // Créer un nom de fichier
        char filename[50];
        sprintf(filename, "Résultats/similarite/similarite_%d.txt", liste[i]->Id);
        // Ouvrir le fichier en écriture
        FILE *fichier = fopen(filename, "w");
        if (fichier == NULL) {
            fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
            exit(EXIT_FAILURE);
        }

        //Enregistrement
        for (int y = 0; y < 100; y++)
        {
            if (tab[y] != 0.0) {// Écriture dans le fichier
            //printf("\n%d\t%lf",y,tab[y]);
            fprintf(fichier, "Graphe moléculaire d'ID chEBI : %d\nDegré de similarité : %lf", sliste[y]->Id, tab[y]);
            fprintf(fichier, "\nTaille de la base de cycle : %d\t",sliste[y]->nb_cycles);
            fprintf(fichier, "\ttaille des cycles de la base : ");
            for(int x = 0; x < sliste[y]->nb_cycles; x++){
                fprintf(fichier, " %d ,",sliste[y]->generateur[x].taille);
            }
            fprintf(fichier, "\n\n");}
        }
        // Fermer le fichier
        fclose(fichier); 
        free(tab); 
        free(sliste);
    }

}

//retrouver les 100 premiers graphes les plus similaire au graphe passé en paramètres
void graphesSimilaires(struct g_cycles g1){
  printf("\n\n--------------ID CHEBI DE LA MOLÉCULE QUERY : %d --------------\n\n",g1.Id);
    char filename[50];
    sprintf(filename, "Résultats/similarite/similarite_%d.txt", g1.Id);
    // Ouvrir le fichier en mode lecture
    FILE *fichier = fopen(filename, "r");
    int TAILLE_MAX = 100;
    char chaine[TAILLE_MAX];
    if (fichier != NULL)
    {
        //Modification pour afficher juste les 100 premiers
        while (fgets(chaine, TAILLE_MAX, fichier) != NULL) // On lit le fichier tant qu'on ne reçoit pas d'erreur (NULL)
        {
            printf("%s", chaine); // On affiche la chaîne qu'on vient de lire
            
        }
        fclose(fichier);
    }
}

int rechercheIndiceAvecId(struct g_cycles **TabGrapheCycle, int cpt, int id){
    for(int i = 3;i<cpt;i++){
        if(TabGrapheCycle[i]->Id == id){
            return i;
        }
    }
    printf("PROBLEME RECHERCHE INDICE !!!\n");
    return 0;
    }

void classeEquivalences(struct g_cycles **TabGrapheCycle, int cpt){// cpt la taille du tableau
    int ** ClasseEquivalence;
    ClasseEquivalence = malloc(cpt*sizeof(int*));
    for(int i = 0 ; i<cpt;i++){
        ClasseEquivalence[i] = malloc(cpt*sizeof(int));
    }
    int * indiceClasse;
    indiceClasse = malloc(cpt*sizeof(int));
    int nbClasseEqui = 0;
    for(int i = 3;i<cpt;i++){

        int ajout = 0;
        for(int j=0;j<nbClasseEqui;j++){
            if(comparaisonTailleDeCycles(*TabGrapheCycle[i],*TabGrapheCycle[rechercheIndiceAvecId(TabGrapheCycle,cpt,ClasseEquivalence[j][0])])==1){
                ClasseEquivalence[j][indiceClasse[j]] = TabGrapheCycle[i]->Id;
                indiceClasse[j]++;
                ajout = 1;
            }
        }
        if(ajout ==0){
            ClasseEquivalence[nbClasseEqui][0] = TabGrapheCycle[i]->Id;
            indiceClasse[nbClasseEqui] = 1;
            nbClasseEqui++;
        }
    }

     char filename[50];
      sprintf(filename, "Résultats/classesEquivalences.txt");
      // Ouvrir le fichier en écriture
      FILE *fichier = fopen(filename, "w");
      if (fichier == NULL) {
          fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
          exit(EXIT_FAILURE);
      }

    for (int i = 0; i < nbClasseEqui; i++) {
        fprintf(fichier, "Classe %d : ", i);
        int k = 3;
        int id1 = ClasseEquivalence[i][0];
        while (TabGrapheCycle[k]->Id != id1 && k < cpt) {
            k += 1;
        }
        fprintf(fichier, "Taille de la base de cycles : %d, Taille des cycles :  ",
            TabGrapheCycle[k]->nb_cycles);
        for (int j = 0; j < TabGrapheCycle[k]->nb_cycles; j++) {
            fprintf(fichier, "%d, ", TabGrapheCycle[k]->generateur[j].taille);
        }
        fprintf(fichier, "\nListe des IDs dans cette classe d'équivalence : ");
        for (int j = 0; j < indiceClasse[i]; j++) {
            fprintf(fichier, "%d ", ClasseEquivalence[i][j]);
        }
        fprintf(fichier, "\n\n");
    }

    // Libération de la mémoire
    for (int i = 0; i < cpt; i++) {
        free(ClasseEquivalence[i]);
    }
    free(ClasseEquivalence);
    free(indiceClasse);
}