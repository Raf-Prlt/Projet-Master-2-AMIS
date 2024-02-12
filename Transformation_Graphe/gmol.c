#include "gmol.h"
#include <assert.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ajouterVoisin(struct Atome *atome, int voisinId) {
  atome->voisinsIds =
      realloc(atome->voisinsIds, (atome->degre + 1) * sizeof(int));
  atome->voisinsIds[atome->degre] = voisinId;
  atome->degre++;
}

struct Atome getAtomeById(struct g_mol *molecule, int nb_atomes, int id) {
  struct Atome atomeVide = {0}; // Une structure avec des valeurs par défaut

  for (int i = 0; i < nb_atomes; i++) {
    struct Atome a = molecule->Atomes[i];
    if (a.Id == id) {
      return a; // Retourner l'atome si son ID correspond
    }
  }

  return atomeVide; // Retourner une structure spéciale pour indiquer l'absence
                    // d'atome
}

// Fonction pour écrire la structure g_mol dans un fichier
void write_g_mol(struct g_mol *molecule) {
#ifdef _WIN32
  system("mkdir g_mol 2> nul"); // Pour Windows
#else
  system("mkdir -p g_mol 2> /dev/null"); // Pour Linux/Unix
#endif

  // Créer un nom de fichier avec le format "g_mol_(id).dat" pour les données
  // binaires
  char filename[50];
  sprintf(filename, "g_mol/g_mol_%d.dat", molecule->Id);

  // Ouvrir le fichier en écriture binaire
  FILE *file = fopen(filename, "wb");

  if (file == NULL) {
    fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
    exit(EXIT_FAILURE);
  }

  // Écriture de la structure g_mol dans le fichier
  fwrite(molecule, sizeof(struct g_mol), 1, file);

  // Fermer le fichier
  fclose(file);
}

// Fonction pour lire la structure g_mol depuis un fichier
struct g_mol read_g_mol(char *filename) {

  FILE *file = fopen(filename, "rb");

  if (file == NULL) {
    fprintf(stderr, "Erreur lors de l'ouverture du fichier %s\n", filename);
    exit(EXIT_FAILURE);
  }

  // Lecture de la structure g_mol depuis le fichier
  struct g_mol test;
  fread(&test, sizeof(struct g_mol), 1, file);

  // Affichage des données lues
  printf("Name: %s , ID: %d\n", test.name, test.Id);

  // Fermer le fichier
  fclose(file);

  return test;
}

void freeMolecule(struct g_mol *molecule, int nb_atomes) {
  for (int i = 0; i < nb_atomes; i++) {
    free(molecule->Atomes[i].voisinsIds);
  }
  free(molecule->Atomes);
  free(molecule->liaisons);
  free(molecule->name);
  free(molecule);
}

struct g_mol *gmol(FILE *file) {

  // Transformation en graphe moléculaire

  // Pour ignorer la première ligne et la ligne contenant "RDKit 2D"
  char buffer[100];
  fgets(buffer, sizeof(buffer), file);
  fgets(buffer, sizeof(buffer), file);

  int nb_atomes, nb_liaisons;
  long int position_before =
      ftell(file); // Enregistrez la position actuelle du curseur

  // Essayer de lire avec la première spécification de format
  fscanf(file, "%3d%3d", &nb_atomes, &nb_liaisons);

  if (nb_atomes == 0 || nb_liaisons == 0) {
    // La lecture n'a pas réussi à obtenir deux nombres ou les nombres sont nuls
    // printf("Problème lecture: Il y a %d atomes, et %d liaisons\n",nb_atomes,
    // nb_liaisons);
    // Remettre le curseur de fichier au début du champ de données
    fseek(file, position_before, SEEK_SET);

    // Essayer de lire avec la deuxième spécification de format
    fgets(buffer, sizeof(buffer), file);
    fgets(buffer, sizeof(buffer), file);
    fscanf(file, " %2d%3d", &nb_atomes, &nb_liaisons);
  }
  fgets(buffer, sizeof(buffer), file);
  /*
  if(!nb_atomes && !nb_liaisons)
  {
      printf("Fichier %18s:Pas d'atome ni de liaison\n", entry->d_name);
  }

  else if (!nb_atomes)
  {
      printf("Fichier %18s:Pas d'atome\n", entry->d_name);
  }
  else if (!nb_liaisons)
  {
      printf("Fichier %18s:Pas de liaison\n", entry->d_name);
  } */

  printf("Il y a %d atomes, et %d liaisons\n", nb_atomes, nb_liaisons);
  struct Atome *atomes = NULL;
  atomes = malloc(nb_atomes * sizeof(struct Atome));
  struct Liaison *liaisons = NULL;
  liaisons = malloc(nb_liaisons * sizeof(struct Liaison));
  int i = 0;
  // Enregistrement des atomes
  while (i < nb_atomes) {
    atomes[i].Id = i + 1;
    fscanf(file, "%*f %*f %*f %2s", atomes[i].symbole);
    // printf("Symbole de l'atome n°%d : %s\n",atomes[i].Id, atomes[i].symbole);
    fgets(buffer, sizeof(buffer), file);

    i++;
  }

  // Enregistrement des liaisons
  int j = 0;

  while (j < nb_liaisons) {
    long int position_before =
        ftell(file); // Pour avoir la position actuelle du curseur

    // Essayer de lire les indices
    fscanf(file, "\n%3d%3d %d", &liaisons[j].IdA1, &liaisons[j].IdA2,
           &liaisons[j].Poids);
    if (liaisons[j].IdA1 == 0 || liaisons[j].IdA2 == 0 ||
        liaisons[j].Poids == 0 || liaisons[j].IdA1 > nb_atomes ||
        liaisons[j].IdA2 > nb_atomes) {
      // La lecture n'a pas réussi ou les indices ne sont pas valides, on
      // revient au début de la ligne
      fseek(file, position_before, SEEK_SET);
      fscanf(file, "\n %2d%3d %d", &liaisons[j].IdA1, &liaisons[j].IdA2,
             &liaisons[j].Poids);
    }
    fgets(buffer, sizeof(buffer), file);

    // printf("Liaison n°%d: entre les atomes n°%3d et n°%3d, de type %d\n",
    // j+1, liaisons[j].IdA1, liaisons[j].IdA2, liaisons[j].Poids);
    j++;
  }

  struct g_mol *molecule = NULL;
  molecule = malloc(sizeof(struct g_mol));
  molecule->nb_atomes = nb_atomes;
  molecule->nb_liaisons = nb_liaisons;
  molecule->Atomes = atomes;
  molecule->liaisons = liaisons;

  // Remplissage des voisins
  for (int i = 0; i < nb_atomes; i++) {
    molecule->Atomes[i].degre = 0;
    molecule->Atomes[i].voisinsIds = NULL;
  }
  // for (int j = 0; j < nb_liaisons; j++) {
  //   molecule->Atomes[molecule->liaisons[j].IdA1 - 1].degre++;
  //   molecule->Atomes[molecule->liaisons[j].IdA2 - 1].degre++;
  // }

  for (int i = 0; i < nb_liaisons; i++) {
    ajouterVoisin(&atomes[liaisons[i].IdA1 - 1], liaisons[i].IdA2);
    ajouterVoisin(&atomes[liaisons[i].IdA2 - 1], liaisons[i].IdA1);
  }

  char *balise_id = "<ChEBI ID>"; // Balise de l'id à chercher
  molecule->Id = -1;              // Valeur par défaut si non trouvée

  while (fgets(buffer, sizeof(buffer), file) != NULL) {
    // Vérifie si la ligne contient la balise
    if (strstr(buffer, balise_id) != NULL) {
      fscanf(file, "%*[^0123456789]%d", &molecule->Id);
      // printf("ChEBI ID: %d\n", molecule->Id);

      break;
    }
  }

  if (molecule->Id == -1) {
    // printf("Erreur: ChEBI ID non trouvée\n");
  }

  char *balise_name = "<ChEBI Name>"; // Balise du nom à chercher
  while (fgets(buffer, sizeof(buffer), file) != NULL) {
    // Vérifie si la ligne contient la balise
    if (strstr(buffer, balise_name) != NULL) {

      char *start = strchr(buffer, '>'); // Trouve le premier '>'
      fgets(buffer, sizeof(buffer), file);
      if (start != NULL) {
        char *end = strchr(start, '\n');
        if (end != NULL) {
          *end = '\0';

          molecule->name = strdup(start);
        }
      }
      if (!molecule->name) {
        // printf("Erreur: ChEBI Name non trouvée\n");
      }
      // printf("ChEBI name: %s\n", molecule->name);

      break;
    }
  }
  return molecule;
}
