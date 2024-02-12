#ifndef GMOL
#define GMOL

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>

struct Atome {
  int Id;          // Numéro dans le fichier sdf
  char symbole[2]; // trouver un moyen de transformer ça en num atomique
  int degre;
  int *voisinsIds;
  // liaisons, autres infos...
};

struct g_mol {
  int Id; // ChEBI ID
  char *name;
  int nb_atomes;
  int nb_liaisons;
  struct Atome *Atomes;
  struct Liaison *liaisons;
};

struct Liaison {
  int IdA1;
  int IdA2;
  int Poids;
};

struct g_mol *gmol(FILE *file);
void freeMolecule(struct g_mol *molecule, int nb_atomes);
void write_g_mol(struct g_mol *molecule);
struct g_mol read_g_mol(char *filename);

#endif