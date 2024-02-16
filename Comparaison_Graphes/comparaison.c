#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <stdbool.h>

typedef struct atome {
  int Id; // Numéro dans le fichier sdf
  int NumeroAtomique;
  char symbole[3]; // trouver un moyen de transformer ça en num atomique
  int degre;
  struct atome *voisins;
  // liaisons, autres infos...
} Atome;

typedef struct Liaison {
  int IdA1;
  int IdA2;
  int Poids;
} Liaison;

typedef struct g_mol {
  int Id; // ChEBI ID
  char *name;
  int nb_atomes;
  int nb_liaisons;
  Atome *Atomes;
  Liaison *Liaisons;
} G_Mol;

typedef struct Cycle {
  int Id;
  int taille; // nombres d'arêtes
  int *liaisons;
  int degre;
  int *voisinsIds;
} Cycle;

typedef struct LiaisonCycles {
  int IdC1;
  int IdC2;
  int Poids;
} LiaisonCycles;

typedef struct g_cycles {
  int Id;
  int nb_cycles;
  int nb_liaisons;
  Cycle *generateur; // classer par ordre croissant
  LiaisonCycles *aretes;
  G_Mol *molecule;
} Graphe_Cycle;

// Fonction pour comparer deux graphes de cycles sous le critère du nombre de cycles
bool comparaisonNbreDeCycles(Graphe_Cycle g1, Graphe_Cycle g2){

    // regarder aussi le nbre de liaisons?
    if(g1.nb_cycles == g2.nb_cycles) {
        return true;
    }

    return false;
}

// Fonction pour comparer deux graphes de cycles sous le critère de la taille des cycles
double comparaisonTailleDeCycles(Graphe_Cycle g1, Graphe_Cycle g2){
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
double comparaisonNbreVoisinsDeCycle(Graphe_Cycle g1, Graphe_Cycle g2){
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
bool comparaisonPoidsDesLiaisinsDeCycles(Graphe_Cycle g1, Graphe_Cycle g2){

    if(comparaisonNbreVoisinsDeCycle(g1, g2)){
        int poids1 = 0;
        int poids2 = 0;
        for (int i = 0; i < g1.nb_liaisons; i++) {
            poids1 += g1.aretes[i].Poids;
            poids2 += g2.aretes[i].Poids;
        }

        poids1 = poids1/g1.nb_cycles;
        poids2 = poids2/g2.nb_liaisons;

        if (poids1 == poids2) {
            return true;
        }       
    }  

    return false;  
}

// Fonction principale
int main() {
    
    /* - Calculer la similarité entre tous les molécules deux par deux et stocker le résultat 
        dans un fichier sous forme de matrice 
       - Une fonction qui prend en paramètres un graphe de cycle et retourne la liste des 100 premiers
        graphes de cycles qui sont le proches en similarité (utiliser les données stockées dans le fichier) 
       - Ajouter les scores de similarité 
       - Installer la librairie  NAUTY
    */
    return 0;
}