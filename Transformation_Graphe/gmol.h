#include <stdio.h>
#include <stdlib.h>

/*
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
*/

/*
Structure sparsegraph:
 As described in Section 3, the sparse representation of a graph uses a structure of type
 sparsegraph with the following fields:
 int nv: the number of vertices
 size_t nde: the number of directed edges (loops count as 1)
 size_t ∗v: pointer to an array of length at least nv
 int ∗d: pointer to an array of length at least nv
 int ∗e: pointer to an array of length at least nde
 sg_weight ∗w: not implemented in this version, should be NULL
 size_t vlen, dlen, elen, wlen: the actual lengths of the arrays v, d, e and w. The unit
 is the element type of the array in each case (so vlen is the number of ints in the
 array v, etc.)

*/

struct Molecular_graph {
    int nv;                 // (nv) Nombre d'atomes 
    int nde;                // (nde) Nombre de liaisons
    int *d;                 // (d) Tableau des degrés pour chaque atome
    int *liaison_id1;       // Tableau des indices du 1er atome pour chaque liaison
    int *liaison_id2;       // Tableau des indices du 2nd atome pour chaque liaison
    int *w;                 // (w) Tableau des types de liaisons pour chaque arête
    int **e;                // (e) Tableau des indices des voisins pour chaque atome
    char **symb_atom;       // Tableau de chaînes de caractères pour les symboles des atomes
    char *chebi_name; 
    int chebi_id;           
};
