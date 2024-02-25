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
        while (j < g2.nb_cycles)
        {
            if (g1.generateur[i].taille == g2.generateur[j].taille){
                cyclesIdentiques++;
                j++;
            }else if (g1.generateur[i].taille < g2.generateur[j].taille)
            {
                if (i == g1.nb_cycles)
                {
                    j = g2.nb_cycles;
                }else{
                    i++;
                }
            }else if (g1.generateur[i].taille > g2.generateur[j].taille)
            {
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
    degreDeSimilarite = (double) min / max;
    return degreDeSimilarite;
}

// Fonction pour comparer deux graphes de cycles sous le critère du poids moyen des liaisons entre les cycles
double comparaisonPoidsDesArretes(struct g_cycles g1, struct g_cycles g2){
        double degreDeSimilarite;
        double poids1 = 0;
        double poids2 = 0;
        double max;
        double min;

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
    degreDeSimilarite = (double) min / max;
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
void calculSimilarite(struct g_cycles *liste, int size){
    //Création du sous dossier
    #ifdef _WIN32
        system("mkdir similarite 2> nul"); // Pour Windows
    #else
        system("mkdir -p g_mol 2> /dev/null"); // Pour Linux/Unix
    #endif

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

    for (int i = 0; i < size; i++)
    {
        //Définition d'une sous liste
        int *sliste = malloc(size * sizeof(int));
        for (int k = 0; k < size; k++)
        {
            sliste[k] = liste[k].Id;
        }

        // Tri décroissant par ligne
        for (int k = 0; k < size-1; k++)
        {
            for (int j = 1; j < size; j++)
            { 
                if(tab[i][k] < tab[i][j])
                {
                    int x = sliste[k];
                    double c = tab[i][k];

                    tab[i][k] = tab[i][j];
                    sliste[k] = sliste[j];

                    tab[i][j] = c;
                    sliste[j] = x;    
                }
            } 
        }

        // Créer un nom de fichier
        char filename[50];
        sprintf(filename, "similarite/similarite_%d.txt", liste[i].Id);
        // Ouvrir le fichier en écriture
        FILE *fichier = fopen(filename, "w");
        if (fichier == NULL) {
            fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
            exit(EXIT_FAILURE);
        }

        //Enregistrement
        for (int j = 0; j < size; j++)
        {
            // Écriture dans le fichier
            //fprintf(fichier, "%d %d %s %lf", liste[j].Id, liste[j].molecule->Id, liste[j].molecule->name, tab[i][j]);
            fprintf(fichier, "%lf  %d", tab[i][j], sliste[j]);
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


// Fonction principale pour faire des tests
int mainComp() {
    /*
    struct g_cycles *liste = malloc(3*sizeof(struct g_cycles));

    struct g_cycles g1;
    struct g_cycles g2;
    struct g_cycles g3;
    
    
    //Test nombre de cycle
    printf("Similarité nombre de cycle: %lf \n", comparaisonNbreDeCycles(g1, g2));
    printf("Similarité nombre de cycle: %lf \n", comparaisonNbreDeCycles(g1, g3));

    //Test taille de cycle
    printf("Similarité taille de cycle: %lf \n", comparaisonTailleDeCycles(g1, g2));
    printf("Similarité taille de cycle: %lf \n", comparaisonTailleDeCycles(g1, g3));

    //Test nombre de voisins
    printf("Similarité nombre de voisins: %lf \n", comparaisonNbreVoisinsDeCycle(g1, g2));
    printf("Similarité nombre de voisins: %lf \n", comparaisonNbreVoisinsDeCycle(g1, g3));

    //Test poids des arêtes
    printf("Similarité poids des arêtes: %lf \n", comparaisonPoidsDesArretes(g1, g2));
    printf("Similarité poids des arêtes: %lf \n", comparaisonPoidsDesArretes(g1, g3));

    //Test similarité
    printf("Degré de similarité : %lf \n", similarite(g1, g2));
    printf("Degré de similarité : %lf \n", similarite(g1, g3));
   
    //Test calcul des similarités
    calculSimilarite(liste, 3);
    graphesSimilaires(g1);
    */
   return 0;
}