#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <stdbool.h>
#include "comparaison.h"

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de cycles
bool comparaisonNbreDeCycles(struct g_cycles g1, struct g_cycles g2){

    // regarder aussi le nbre de liaisons?
    if(g1.nb_cycles == g2.nb_cycles) {
        return true;
    }

    return false;
}

// Fonction pour comparer deux graphes de cycles sous le critère de la taille des cycles
double comparaisonTailleDeCycles(struct g_cycles g1, struct g_cycles g2){
    int cyclesIdentiques = 0;
    if(comparaisonNbreDeCycles(g1, g2)){
        for (int i = 0; i < g1.nb_cycles; i++)
        {
            if(g1.generateur[i].taille != g2.generateur[i].taille) {
                i = g1.nb_cycles;
            }else{
                cyclesIdentiques++;
            }
        }
        return (1 + (cyclesIdentiques/g1.nb_cycles));
    } else {
        int max;
        int min;
        if (g1.nb_cycles < g2.nb_cycles){
            min = g1.nb_cycles;
            max = g2.nb_cycles;
        }else{
            max = g1.nb_cycles;
            min = g2.nb_cycles;
        }

        for (int i = 0; i < min; i++){
            if (g1.generateur[i].taille != g2.generateur[i].taille){
                i = min;
            }else{
                cyclesIdentiques++;
            } 
        }   
        return cyclesIdentiques/max;     
        
    }
    return 0;
}

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de voisins
double comparaisonNbreVoisinsDeCycle(struct g_cycles g1,struct g_cycles g2){
    int cyclesDegreIdentiques = 0;
    if (comparaisonTailleDeCycles(g1, g2) == 2){
        for (int i = 0; i < g1.nb_cycles; i++) {
            if (g1.generateur[i].degre == g2.generateur[i].degre)
            {
                cyclesDegreIdentiques++;
            }
            
        }
        return  2+(cyclesDegreIdentiques/g1.nb_cycles);
    }else{
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

        if (degre1 == degre2) {
            return 0.8 + comparaisonTailleDeCycles(g1, g2);
        }else if ((degre1 < (degre2+1)) && ((degre2-1) <= degre1)){
            return 0.5 + comparaisonTailleDeCycles(g1, g2);
        }else if ((degre1 < (degre2+2)) && ((degre2-2) <= degre1)){
            return 0.1 + comparaisonTailleDeCycles(g1, g2);
        }else{
            return comparaisonTailleDeCycles(g1, g2);
        }
    }
    return 0;  
}

// Fonction pour comparer deux graphes de cycles sous le critère du poids moyen des liaisons entre les cycles
double degreDeSimilarite(struct g_cycles g1, struct g_cycles g2){

        int poids1 = 0;
        int poids2 = 0;
        for (int i = 0; i < g1.nb_liaisons; i++) {
            poids1 += g1.aretes[i].Poids;
            poids2 += g2.aretes[i].Poids;
        }

        poids1 = poids1/g1.nb_cycles;
        poids2 = poids2/g2.nb_liaisons;

        if (poids1 == poids2) {
            return 1 + comparaisonNbreVoisinsDeCycle(g1, g2);
        }else if ((poids1 < (poids2+1)) && ((poids2-1) <= poids1)){
            return 0.5 + comparaisonNbreVoisinsDeCycle(g1, g2);
        }else if ((poids1 < (poids2+2)) && ((poids2-2) <= poids1)){
            return 0.1 + comparaisonNbreVoisinsDeCycle(g1, g2);
        }else{
            return comparaisonNbreVoisinsDeCycle(g1, g2);
        }    

    return 0;
}

// Calculer les similarités entre les graphes de cycles et stocker le résultat dans un fichier
void calculSimilarites(struct g_cycles *liste){
    double tab[5][5];

    for (int i = 0; i < 5; i++)
    {
        for (int j = i; j < 5; j++)
        {
            tab[i][j] = degreDeSimilarite(liste[i], liste[j]);
            tab[j][i] = tab[i][j];
        }
        
    }

    FILE* fichier = NULL;
    fichier = fopen("fichier.txt", "r+");

    if (fichier != NULL)
    {
        for (int i = 0; i < 5; i++){
            fprintf(fichier, "%d ", liste[i].Id);
            for (int j = 0; j < 5; j++)
            {
                fprintf(fichier, "%d ", tab[i][j]);
            }
            fprintf(fichier, "\n");
        }
    }

    fclose(fichier);
    
}

//retrouver les 100 premiers graphes les plus similaire au graphe passé en paramètres
void similarites(struct g_cycles *liste, struct g_cycles g1){



}


// Fonction principale pour faire des tests
int mainComp() {
    int tab[5][5];

    for (int i = 0; i < 5; i++)
    {
        for (int j = i; j < 5; j++)
        {
            tab[i][j] = i;
            tab[j][i] = i;
        }
        
    }

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            printf("%d ", tab[i][j]);
        }
        printf("\n");
        
    }

    FILE* fichier = NULL;
    fichier = fopen("fichier.txt", "r+");

    if (fichier != NULL)
    {
        for (int i = 0; i < 5; i++){
            for (int j = 0; j < 5; j++)
            {
                //  fputc(tab[i][j]+" ", fichier);
                fprintf(fichier, "%d ", tab[i][j]);
            }
            fprintf(fichier, "\n");
        }
        printf("OK \n");
    }

    fclose(fichier);
    
    
    /* - Calculer la similarité entre tous les molécules deux par deux et stocker le résultat 
        dans un fichier sous forme de matrice 
       - Une fonction qui prend en paramètres un graphe de cycle et retourne la liste des 100 premiers
        graphes de cycles qui sont le proches en similarité (utiliser les données stockées dans le fichier) 
       - Ajouter les scores de similarité 
       - Installer la librairie  NAUTY
    */
    return 0;
}