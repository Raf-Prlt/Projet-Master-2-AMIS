#include <stdio.h>
#include <stdlib.h>

struct Atome{
    int Id;             //Numéro dans le fichier sdf
    int NumeroAtomique;
    char symbole[3];    // trouver un moyen de transformer ça en num atomique 
    struct Atome* voisins;    
    // liaisons, autres infos...
};

struct g_mol{
    int Id;         //ChEBI ID
    char* name;    
    struct Atome* Atomes; 
};

struct Liaison{
    int IdA1;   
    int IdA2;   
    int Poids;  
};