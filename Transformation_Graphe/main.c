#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "gmol.h"

int main() {
    DIR *directory;
    struct dirent *entry;
    directory = opendir("../TraitementDonnee/test/"); // à modifier quand ce sera bon (remplacer test par output_molecules)
    if (!directory) {
        printf("Impossible d'ouvrir le dossier");
        exit(1);
    }

    // Lit chaque fichier dans le dossier
    while ((entry = readdir(directory)) != NULL) {

        // Vérifie si le fichier est un fichier .sdf
        if (strstr(entry->d_name, ".sdf") != NULL) {

            // Construit le chemin complet du fichier
            char filepath[512];
            snprintf(filepath, sizeof(filepath), "../TraitementDonnee/test/%s", entry->d_name); // à modifier quand ce sera bon
            FILE *file = fopen(filepath, "r");
            if (!file) {
                printf("Impossible d'ouvrir le fichier");
            }
            
            //Transformation en graphe moléculaire

            // Pour ignorer la première ligne et la ligne contenant "RDKit 2D"
            char buffer[100];
            fgets(buffer, sizeof(buffer), file);
            fgets(buffer, sizeof(buffer), file);
            
            int nb_atomes, nb_liaisons;
            fscanf(file, "%d %d", &nb_atomes, &nb_liaisons);
            
            fgets(buffer, sizeof(buffer), file);

            printf("Il y a %d atomes, et %d liaisons\n",nb_atomes, nb_liaisons);
            struct Atome *atomes = malloc(nb_atomes * sizeof(struct Atome));
            struct Liaison *liaisons = malloc(nb_liaisons * sizeof(struct Liaison));
            int i = 0;
            // Enregistrement des atomes
            while (i < nb_atomes) {
                atomes[i].Id = i+1;
                fscanf(file, "%*f %*f %*f %s", atomes[i].symbole);
                printf("Symbole de l'atome n°%d : %s\n",atomes[i].Id, atomes[i].symbole);
                fgets(buffer, sizeof(buffer), file);

                i++;
            }

            // Enregistrement des liaisons
            int j = 0;
            while (j < nb_liaisons) {
                fscanf(file, "%d %d %d", &liaisons[j].IdA1, &liaisons[j].IdA2, &liaisons[j].Poids);
                printf("Liaison entre les atomes n°%2d et n°%2d, de type %d\n",liaisons[j].IdA1, liaisons[j].IdA2, liaisons[j].Poids);
                fgets(buffer, sizeof(buffer), file);

                j++;
            }

            struct g_mol* molecule = malloc(sizeof(struct g_mol));
            molecule->Atomes = atomes;

            char* balise_id = "<ChEBI ID>"; // Balise de l'id à chercher
            molecule->Id = -1; // Valeur par défaut si non trouvée

            while (fgets(buffer, sizeof(buffer), file) != NULL) {
                // Vérifie si la ligne contient la balise
                if (strstr(buffer, balise_id) != NULL) {
                    fscanf(file, "%*[^0123456789]%d", &molecule->Id);
                    printf("ChEBI ID: %d\n", molecule->Id);

                    break;
                }
            }

            if (molecule->Id == -1) {
                printf("Erreur: Balise non trouvée\n");
            }
            
            char* balise_name = "<ChEBI Name>"; // Balise du nom à chercher
            while (fgets(buffer, sizeof(buffer), file) != NULL) {
                // Vérifie si la ligne contient la balise
                if (strstr(buffer, balise_name) != NULL) {

                    char* start = strchr(buffer, '>'); // Trouve le premier '>'
                    fgets(buffer, sizeof(buffer), file);
                    if (start != NULL) {
                        char* end = strchr(start, '\n'); 
                        if (end != NULL) {
                            *end = '\0';

                            molecule->name = strdup(start);
                        }
                    }
                    printf("ChEBI name: %s\n", molecule->name);

                    break;
                }
            }

            
            free(atomes);
            free(liaisons);
            free(molecule->name);
            free(molecule);
            fclose(file);
        }
    }

    closedir(directory);
    return 0;
}
