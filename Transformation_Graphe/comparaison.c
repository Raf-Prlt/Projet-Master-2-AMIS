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
    degreDeSimilarite = min/max;

    return degreDeSimilarite;
}

// Fonction pour comparer deux graphes de cycles sous le critère de la taille des cycles
double comparaisonTailleDeCycles(struct g_cycles g1, struct g_cycles g2){
        double degreDeSimilarite;
        int cyclesIdentiques = 0;
        int max;
        int min;
        if (g1.nb_cycles < g2.nb_cycles){
            min = g1.nb_cycles;
            max = g2.nb_cycles;
        }else{
            max = g1.nb_cycles;
            min = g2.nb_cycles;
        }
        int i = 0;
        int j = 0;
        while (i < g1.nb_cycles)
        {
            while (j < g2.nb_cycles)
            {
                if (g1.generateur[i].taille == g2.generateur[i].taille){
                    cyclesIdentiques++;
                }else if (g1.generateur[i].taille < g2.generateur[i].taille)
                {
                    i++;
                }else if (g1.generateur[i].taille > g2.generateur[i].taille)
                {
                    j++;
                }  
            }
        }
        degreDeSimilarite = cyclesIdentiques/max;
        return degreDeSimilarite;  
}

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de voisins
double comparaisonNbreVoisinsDeCycle(struct g_cycles g1,struct g_cycles g2){
    double degreDeSimilarite;
    int max;
    int min;
    int degre1 = 0;
    int degre2 = 0;

    for (int i = 0; i < g1.nb_cycles; i++) {
        degre1 += g1.generateur[i].degre;
    }
    for (int i = 0; i < g2.nb_cycles; i++) {
        degre2 += g2.generateur[i].degre;
    }

    degre1 = degre1/g1.nb_cycles;
    degre2 = degre2/g2.nb_cycles;

    if (degre1 < degre2){
        min = degre1;
        max = degre2;
    }else{
        max = degre1;
        min = degre2;
    }
    degreDeSimilarite = min/max;
    return degreDeSimilarite;
}

// Fonction pour comparer deux graphes de cycles sous le critère du poids moyen des liaisons entre les cycles
double comparaisonPoidsDesArretes(struct g_cycles g1, struct g_cycles g2){
        double degreDeSimilarite;
        int poids1 = 0;
        int poids2 = 0;
        int max;
        int min;

        for (int i = 0; i < g1.nb_liaisons; i++) {
            poids1 += g1.aretes[i].Poids;
            poids2 += g2.aretes[i].Poids;
        }

        poids1 = poids1/g1.nb_cycles;
        poids2 = poids2/g2.nb_liaisons; 

        if (poids1 < poids2){
            min = poids1;
            max = poids2;
        }else{
            max = poids1;
            min = poids2;
        }
    degreDeSimilarite = min/max;
    return degreDeSimilarite;
}

double similarite(struct g_cycles g1, struct g_cycles g2){
    double degreDeSimilarite;
    double degre = 0;

    degre += comparaisonNbreDeCycles(g1, g2);
    degre += comparaisonTailleDeCycles(g1, g2);
    degre += comparaisonNbreVoisinsDeCycle(g1, g2);
    degre += comparaisonPoidsDesArretes(g1, g2);

    degreDeSimilarite = degre/4;
    return degreDeSimilarite;

}

// Calculer les similarités entre les graphes de cycles et stocker le résultat dans un fichier
void calculSimilarite(struct g_cycles *liste, int size){
    double **tab = malloc(size * sizeof(double*));
    
    for (int i = 0; i < size; i++)
    {
        tab[i] = malloc(size * sizeof(double));
    }
    // Calcul du degré de similarité
    for (int i = 0; i < size; i++)
    {
        tab[i][i] = 0; // pas de comparaison avec soi-même
        for (int j = i+1; j < size; j++)
        {
            tab[i][j] = similarite(liste[i], liste[j]);
            tab[j][i] = tab[i][j];
        } 
    }

    // Enregistrement
    for (int i = 0; i < size; i++)
    {
        // Créer un nom de fichier
        char filename[50];
        sprintf(filename, "similarite/similarite_%d.txt", liste[i].Id);
        // Ouvrir le fichier en écriture
        FILE *fichier = fopen(filename, "w");
        if (fichier == NULL) {
            fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
            exit(EXIT_FAILURE);
        }
        for (int j = 0; j < size; j++)
        {
            // Écriture dans le fichier
            fprintf(fichier, "%d %d %s %lf", liste[j].Id, liste[j].molecule->Id, liste[j].molecule->name, tab[i][j]);
            fprintf(fichier, "\n");
        }
        // Fermer le fichier
        fclose(fichier); 
    }    
}

//retrouver les 100 premiers graphes les plus similaire au graphe passé en paramètres
void graphesSimilaires(struct g_cycles g1){
    char filename[50];
    sprintf(filename, "similarite/similarite_%d.txt", g1.Id);
    // Ouvrir le fichier en mode lecture
    FILE *fichier = fopen(filename, "r");
    int TAILLE_MAX = 1000;
    char chaine[TAILLE_MAX];
    if (fichier != NULL)
    {
        while (fgets(chaine, TAILLE_MAX, fichier) != NULL) // On lit le fichier tant qu'on ne reçoit pas d'erreur (NULL)
        {
            printf("%s", chaine); // On affiche la chaîne qu'on vient de lire
        }
        // A suivre ...
        fclose(fichier);
    }
}


// Fonction principale pour faire des tests
int mainComp() {
    struct g_cycles g1;
    struct g_cycles g2;
    struct g_cycles g3;

    //initialisation
    //Id
    g1.Id = 1;
    g2.Id = 2;
    g3.Id = 3;

    //Nombre de cycles;
    g1.nb_cycles = 3;
    g2.nb_cycles = 3;
    g3.nb_cycles = 4;

    //Test nombre de cycle
    printf("Similarité nombre de cycle: %lf \n", comparaisonNbreDeCycles(g1, g2));
    printf("Similarité nombre de cycle: %lf \n", comparaisonNbreDeCycles(g1, g3));

    return 0;
}