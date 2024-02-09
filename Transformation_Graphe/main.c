#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "gmol.h"

// Fonction pour allouer la mémoire pour la structure
struct Molecular_graph* allocateMolecularGraph(int nv, int nde) {
    struct Molecular_graph* mol_graph = malloc(sizeof(struct Molecular_graph));

    if (mol_graph != NULL) {
        mol_graph->nv = nv;
        mol_graph->nde = nde;

        mol_graph->d = malloc(nv * sizeof(int));
        mol_graph->liaison_id1 = malloc(nde * sizeof(int));
        mol_graph->liaison_id2 = malloc(nde * sizeof(int));
        mol_graph->w = malloc(nde * sizeof(int));
        mol_graph->symb_atom = malloc(nv * sizeof(char*));
        for (int i = 0; i< nv; i++) {
            mol_graph->symb_atom[i] = malloc(2*sizeof(char));
        }

        mol_graph->chebi_name = NULL;  // Pas besoin d'allouer de la mémoire ici
    }

    return mol_graph;
}

// Fonction pour allouer la mémoire pour le champ e
void allocatee(struct Molecular_graph* mol_graph) {
    mol_graph->e = malloc(mol_graph->nv * sizeof(int*));

    for (int i = 0; i < mol_graph->nv; i++) {
        mol_graph->e[i] = malloc(mol_graph->d[i] * sizeof(int));
    }
}


// Fonction pour libérer la mémoire allouée pour la structure
void freeMolecularGraph(struct Molecular_graph* mol_graph) {
    if (mol_graph != NULL) {

        if (mol_graph->e != NULL) {
            for (int i = 0; i < mol_graph->nv; i++) {
                free(mol_graph->e[i]);
            }
            free(mol_graph->e);
        }        
        free(mol_graph->d);
        free(mol_graph->liaison_id1);
        free(mol_graph->liaison_id2);
        free(mol_graph->w);
        for (int i = 0; i < mol_graph->nv; i++) {
            free(mol_graph->symb_atom[i]);
        }
        free(mol_graph->symb_atom);
        free(mol_graph);
    }
}

void printMolecularGraph(const struct Molecular_graph* mol_graph) {
    if (mol_graph != NULL) {
        printf("Nombre d'atomes : %d\n", mol_graph->nv);
        printf("Nombre de liaisons : %d\n", mol_graph->nde);

        printf("Degrés des atomes : ");
        for (int i = 0; i < mol_graph->nv; i++) {
            printf("%d ", mol_graph->d[i]);
        }
        printf("\n");

        printf("Indices des atomes pour chaque liaison : ");
        for (int i = 0; i < mol_graph->nde; i++) {
            printf("(%d, %d) ", mol_graph->liaison_id1[i], mol_graph->liaison_id2[i]);
        }
        printf("\n");

        printf("Types de liaisons : ");
        for (int i = 0; i < mol_graph->nde; i++) {
            printf("%d ", mol_graph->w[i]);
        }
        printf("\n");

        printf("Symboles d'atomes : ");
        for (int i = 0; i < mol_graph->nv; i++) {
            printf("%s ", mol_graph->symb_atom[i]);
        }
        printf("\n");
        /*
        int nb_e = 0;
        printf("Indices des e pour chaque atome : \n");
        for (int i = 0; i < mol_graph->nv; i++) {
            printf("%3d: ",i+1);
            for (int j = 0; j < mol_graph->d[i]; j++){
                printf("%3d ", mol_graph->e[i][j]);
                nb_e++;
            }
            printf("\n");
        }
        printf("\n");*/

        printf("Nom ChEBI : %s\n", mol_graph->chebi_name);
        printf("ID ChEBI : %d\n Vérif e: ", mol_graph->chebi_id);
        /*if(nb_e == mol_graph->nde*2) {printf("Cool !");}
        else {printf("Pas cool");}*/
        printf("\n");
    }
}

int main() {
    DIR *directory;
    struct dirent *entry;
    char *dir_name = "../TraitementDonnee/output_molecules/";
    directory = opendir(dir_name);
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
            snprintf(filepath, sizeof(filepath), "../TraitementDonnee/output_molecules/%s", entry->d_name);
            FILE *file = fopen(filepath, "r");
            if (!file) {
                printf("Impossible d'ouvrir le fichier");
            }

            //Transformation en graphe moléculaire


            // Pour ignorer la première ligne et la ligne contenant "RDKit 2D"
            char buffer[100];
            fgets(buffer, sizeof(buffer), file);
            fgets(buffer, sizeof(buffer), file);

            int nv, nde;
            long int position_before = ftell(file);  // Enregistrez la position actuelle du curseur

            // Essayer de lire avec la première spécification de format
            fscanf(file, "%3d%3d", &nv, &nde);

            if (nv == 0 || nde == 0) {
                // La lecture n'a pas réussi à obtenir deux nombres ou les nombres sont nuls
                //printf("Problème lecture: Il y a %d atomes, et %d liaisons\n",nv, nde);
                // Remettre le curseur de fichier au début du champ de données
                fseek(file, position_before, SEEK_SET);

                // Essayer de lire avec la deuxième spécification de format
                fgets(buffer, sizeof(buffer), file);
                fgets(buffer, sizeof(buffer), file);
                fscanf(file, " %2d%3d", &nv, &nde);
            }
            fgets(buffer, sizeof(buffer), file);
            /*
            if(!nv && !nde)
            {
                printf("Fichier %18s:Pas d'atome ni de liaison\n", entry->d_name);
            }
            
            else if (!nv)
            {
                printf("Fichier %18s:Pas d'atome\n", entry->d_name);
            }
            else if (!nde)
            {
                printf("Fichier %18s:Pas de liaison\n", entry->d_name);
            } */     

            printf("Il y a %d atomes, et %d liaisons\n",nv, nde);


            struct Molecular_graph *graph = allocateMolecularGraph(nv, nde);
            int i = 0;
            // Enregistrement des atomes
            while (i < nv) {
                //atomes[i].Id = i+1;
                fscanf(file, "%*f %*f %*f %s", graph->symb_atom[i]);
                printf("Symbole de l'atome n°%d : %s\n",i+1, graph->symb_atom[i]);
                fgets(buffer, sizeof(buffer), file);

                i++;
            }

            // Enregistrement des liaisons
            int j = 0;
            while (j < nde) {
                long int position_before = ftell(file);  // Pour avoir la position actuelle du curseur

                // Essayer de lire les indices
                fscanf(file, "\n%3d%3d %d", &graph->liaison_id1[j], &graph->liaison_id2[j], &graph->w[j]);
                if(graph->liaison_id1[j] == 0 || graph->liaison_id2[j] == 0 || graph->w[j] == 0
                    || graph->liaison_id1[j] > nv || graph->liaison_id2[j] > nv) {
                    // La lecture n'a pas réussi ou les indices ne sont pas valides, on revient au début de la ligne
                    fseek(file, position_before, SEEK_SET);
                    fscanf(file, "\n %2d%3d %d", &graph->liaison_id1[j], &graph->liaison_id2[j], &graph->w[j]);
                }
                fgets(buffer, sizeof(buffer), file);

                printf("Liaison n°%d: entre les atomes n°%3d et n°%3d, de type %d\n", j+1, graph->liaison_id1[j], graph->liaison_id2[j], graph->w[j]);
                graph->d[graph->liaison_id1[j]-1]++;
                graph->d[graph->liaison_id2[j]-1]++;
                j++;

            }

            allocatee(graph);
            // remplissage du tableau e
            for (int i = 0; i < nv; i++) {
                for (int j = 0; j < graph->d[i]; j++) {
                    graph->e[i][j] = 0;
                }
            }
            int * degre_temp = malloc(nv * sizeof(int));
            for (int i = 0; i < nde; i++) {
                degre_temp[graph->liaison_id1[i]]++;
                degre_temp[graph->liaison_id2[i]]++;
                graph->e[graph->liaison_id1[i]-1][degre_temp[graph->liaison_id1[i]]-1] = graph->liaison_id2[i];
                graph->e[graph->liaison_id2[i]-1][degre_temp[graph->liaison_id2[i]]-1] = graph->liaison_id1[i];
            }
            free(degre_temp);

            char* balise_id = "<ChEBI ID>"; // Balise de l'id à chercher
            graph->chebi_id = -1; // Valeur par défaut si non trouvée

            while (fgets(buffer, sizeof(buffer), file) != NULL) {
                // Vérifie si la ligne contient la balise
                if (strstr(buffer, balise_id) != NULL) {
                    fscanf(file, "%*[^0123456789]%d", &graph->chebi_id);
                    printf("ChEBI ID: %d\n", graph->chebi_id);

                    break;
                }
            }

            if (graph->chebi_id == -1) {
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

                            graph->chebi_name = strdup(start);
                        }
                    }
                    printf("ChEBI name: %s\n", graph->chebi_name);

                    break;
                }
            }


            printMolecularGraph(graph);
            printf("Fermeture fichier %18s\n\n", entry->d_name);
            //freeMolecularGraph(graph);
            fclose(file);
        }

    }

    closedir(directory);
    return 0;
}
