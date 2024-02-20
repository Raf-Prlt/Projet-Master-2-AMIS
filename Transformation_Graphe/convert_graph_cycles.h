#ifndef CONVERT_GRAPH_CYCLES
#define CONVERT_GRAPH_CYCLES

#include "gmol.h"
#include "nauty2_8_8/nausparse.h" /* which includes nauty.h */
#include <stdio.h>
#include <stdlib.h>

typedef struct Chemin {
    int taille;
    int *atomesIds;
    int *liaisons;
} chemin;

struct Cycle {
    int Id;
    int taille; // nombres d'arêtes
    int *liaisons;
    int degre;
    int *voisinsIds;
};

struct LiaisonCycles {
    int IdC1;
    int IdC2;
    int Poids;
};

struct g_cycles {
    int Id;
    int nb_cycles;
    int nb_liaisons;
    struct Cycle *generateur; // nœuds du graphe de cycle
    struct LiaisonCycles *aretes;
    struct g_mol *molecule;
};

void mcKay(struct g_mol *g, sparsegraph *cg);
void print_sg(sparsegraph *sg);
void print_path(chemin *c, int m);
void init_path(chemin *c, int n, int m);
void free_path(chemin *c);
void smallest_paths(sparsegraph *sg, chemin **ps_cts_chms);
void gmol_to_sparsegraph(struct g_mol *g, sparsegraph *sg);
void Horton(sparsegraph *sg, chemin **tab);
void init_liaisons(sparsegraph *sg, struct Liaison A[]);
int verifIntersectionVide(chemin *c1, chemin *c2, int m, int Som);
void ajoutCycles(struct Cycle *Ci, chemin *c1, chemin *c2,struct Liaison a, int compteur, int m, sparsegraph *sg);
void TriCroissant(struct Cycle Ci[], int tailleTab);
int ExtractionBase(struct Cycle *Ci, int m, struct Cycle *Base,  int tailleCi);
void TriCroissantArete(struct Liaison A[], int tailleTab);



#endif