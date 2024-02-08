#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "gmol.h"

// Fonction pour allouer la mémoire pour la structure
struct Molecular_graph* allocateMolecularGraph(int nb_atomes, int nb_liaisons) {
    struct Molecular_graph* mol_graph = malloc(sizeof(struct Molecular_graph));

    if (mol_graph != NULL) {
        mol_graph->nb_atomes = nb_atomes;
        mol_graph->nb_liaisons = nb_liaisons;

        mol_graph->degres = malloc(nb_atomes * sizeof(int));
        mol_graph->liaison_id1 = malloc(nb_liaisons * sizeof(int));
        mol_graph->liaison_id2 = malloc(nb_liaisons * sizeof(int));
        mol_graph->type_laison = malloc(nb_liaisons * sizeof(int));
        mol_graph->symb_atom = malloc(nb_atomes * sizeof(char*));
        mol_graph->chebi_name = NULL;  // Pas besoin d'allouer de la mémoire ici
    }

    return mol_graph;
}

// Fonction pour allouer la mémoire pour le champ voisins
void allocateVoisins(struct Molecular_graph* mol_graph) {
    mol_graph->voisins = malloc(mol_graph->nb_atomes * sizeof(int*));

    for (int i = 0; i < mol_graph->nb_atomes; i++) {
        mol_graph->voisins[i] = malloc(mol_graph->degres[i] * sizeof(int));
    }
}


// Fonction pour libérer la mémoire allouée pour la structure
void freeMolecularGraph(struct Molecular_graph* mol_graph) {
    if (mol_graph != NULL) {

        if (mol_graph->voisins != NULL) {
            for (int i = 0; i < mol_graph->nb_atomes; i++) {
                free(mol_graph->voisins[i]);
            }
            free(mol_graph->voisins);
        }        
        free(mol_graph->degres);
        free(mol_graph->liaison_id1);
        free(mol_graph->liaison_id2);
        free(mol_graph->type_laison);
        free(mol_graph->symb_atom);
        free(mol_graph->chebi_name);
        free(mol_graph);
    }
}

void printMolecularGraph(const struct Molecular_graph* mol_graph) {
    if (mol_graph != NULL) {
        printf("Nombre d'atomes : %d\n", mol_graph->nb_atomes);
        printf("Nombre de liaisons : %d\n", mol_graph->nb_liaisons);

        printf("Degrés des atomes : ");
        for (int i = 0; i < mol_graph->nb_atomes; i++) {
            printf("%d ", mol_graph->degres[i]);
        }
        printf("\n");

        printf("Indices des atomes pour chaque liaison : ");
        for (int i = 0; i < mol_graph->nb_liaisons; i++) {
            printf("(%d, %d) ", mol_graph->liaison_id1[i], mol_graph->liaison_id2[i]);
        }
        printf("\n");

        printf("Types de liaisons : ");
        for (int i = 0; i < mol_graph->nb_liaisons; i++) {
            printf("%d ", mol_graph->type_laison[i]);
        }
        printf("\n");

        printf("Indices des voisins pour chaque atome : \n");
        for (int i = 0; i < mol_graph->nb_atomes; i++) {
            printf("%d: ",i+1);
            for (int j = 0; j < mol_graph->degres[i]; j++){
                printf("%d ", mol_graph->voisins[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        printf("Symboles d'atomes : ");
        for (int i = 0; i < mol_graph->nb_atomes; i++) {
            printf("%c ", mol_graph->symb_atom[i]);
        }
        printf("\n");

        printf("Nom ChEBI : %s\n", mol_graph->chebi_name);
        printf("ID ChEBI : %d\n", mol_graph->chebi_id);
    }
}

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
            // Essayer de lire avec la première spécification de format
            fscanf(file, "%3d%3d", &nb_atomes, &nb_liaisons);

            if (nb_atomes == 0 || nb_liaisons == 0) {
                // La lecture n'a pas réussi à obtenir deux nombres ou les nombres sont nuls
                //printf("Problème lecture\n");
                // Remettre le curseur de fichier au début du champ de données
                fseek(file, -6L, SEEK_CUR);

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
            //printf("Il y a %d atomes, et %d liaisons\n",nb_atomes, nb_liaisons);
            
            struct Molecular_graph *graph = allocateMolecularGraph(nb_atomes, nb_liaisons);

            
            int i = 0;
            // Enregistrement des atomes
            while (i < graph->nb_atomes) {
                //atomes[i].Id = i+1;
                fscanf(file, "%*f %*f %*f %s", &graph->symb_atom[i]);
                //printf("Symbole de l'atome n°%d : %s\n",atomes[i].Id, atomes[i].symbole);
                fgets(buffer, sizeof(buffer), file);

                i++;
            }
            for (int k = 0; k < nb_atomes; k++) {
                graph->degres[k] = 0;
            }

            // Enregistrement des liaisons
            int j = 0;
            while (j < graph->nb_liaisons) {
                int id1, id2, type;
                fscanf(file, "%d %d %d", &id1, &id2, &type);
                //printf("Liaison entre les atomes n°%2d et n°%2d, de type %d\n",liaisons[j].IdA1, liaisons[j].IdA2, liaisons[j].Poids);
                graph->liaison_id1[j] = id1;
                graph->liaison_id2[j] = id2;
                graph->type_laison[j] = type;

                graph->degres[id1-1] ++;
                graph->degres[id2-1] ++;

                fgets(buffer, sizeof(buffer), file);

                j++;
            }

            allocateVoisins(graph);
            // remplissage du tableau voisins
            for (int i = 0; i < nb_atomes; i++) {
                for (int j = 0; j < graph->degres[i]; j++) {
                    graph->voisins[i][j] = 0;
                }
            }
            int * degre_temp = malloc(nb_atomes * sizeof(int));
            for (int i = 0; i < nb_liaisons; i++) {
                degre_temp[graph->liaison_id1[i]]++;
                degre_temp[graph->liaison_id2[i]]++;
                graph->voisins[graph->liaison_id1[i]-1][degre_temp[graph->liaison_id1[i]]-1] = graph->liaison_id2[i];
                graph->voisins[graph->liaison_id2[i]-1][degre_temp[graph->liaison_id2[i]]-1] = graph->liaison_id1[i];
            }
            free(degre_temp);
            
            char* balise_id = "<ChEBI ID>"; // Balise de l'id à chercher
            graph->chebi_id = -1; // Valeur par défaut si non trouvée

            while (fgets(buffer, sizeof(buffer), file) != NULL) {
                // Vérifie si la ligne contient la balise
                if (strstr(buffer, balise_id) != NULL) {
                    fscanf(file, "%*[^0123456789]%d", &graph->chebi_id);
                    //printf("ChEBI ID: %d\n", molecular_graph->chebi_id);

                    break;
                }
            }
            /*
            if (graph->chebi_id == -1) {
                printf("Fichier %18s: Balise ID non trouvée\n", entry->d_name);
            }*/
            
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

                            graph->chebi_name = strdup(start);
                        }
                    }
                    //printf("ChEBI name: %s\n", molecular_graph->chebi_name);

                    break;
                }
            }
            /*
            if(!graph->chebi_name){
                printf("Fichier %18s: Balise name non trouvée\n", entry->d_name);
            }*/
            printMolecularGraph(graph);
            freeMolecularGraph(graph);
        }
    }

    closedir(directory);
    return 0;
}
