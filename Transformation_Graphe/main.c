#include "convert_graph_cycles.h"
#include "gmol.h"
#include "comparaison.h"
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


int main() {
    DIR *directory2;
    directory2 = opendir(
        "../TraitementDonnee/output_molecules_v3000/"); // à modifier quand ce sera
                                                // bon (remplacer test par
                                                // output_molecules)
    if (!directory2) {
        printf("Impossible d'ouvrir le dossier");
        exit(1);
    }
    // Lit chaque fichier dans le dossier

    struct dirent *entry2;
    int cpt = 0;

    while ((entry2 = readdir(directory2)) != NULL) {
      cpt++;
    }

    printf("\nNB FICHIER VU : %d\n",cpt);


    struct g_cycles **TabGrapheCycle;

    TabGrapheCycle = malloc (cpt * sizeof(struct g_cycles));


/*
 ____________________________________________________________________________________
|                                                                                    |
|             LECTURE  ET   CONVERSION   EN   GRAPHE   DE   CYCLES                   |
|____________________________________________________________________________________|

*/


     DIR *directory;
    struct dirent *entry;
    directory = opendir(
        "../TraitementDonnee/output_molecules_v3000/"); // à modifier quand ce sera
                                                // bon (remplacer test par
                                                // output_molecules)
    if (!directory) {
        printf("Impossible d'ouvrir le dossier");
        exit(1);
    }
    int nb_mol_vues = 0;

    while ((entry = readdir(directory)) != NULL && nb_mol_vues<100) {
       /* printf("--------------------------------- Molécule n°%d "
                "---------------------------------\n",
                nb_mol_vues);*/
        nb_mol_vues++;

        // Vérifie si le fichier est un fichier .sdf
        if (strstr(entry->d_name, ".sdf") != NULL) {

            // Construit le chemin complet du fichier
            char filepath[512];
            snprintf(filepath, sizeof(filepath),
                    "../TraitementDonnee/output_molecules_v3000/%s",
                    entry->d_name); // à modifier quand ce sera bon
            FILE *file = fopen(filepath, "r");
            if (!file) {
                printf("Impossible d'ouvrir le fichier");
            }

            // if(strcmp(entry->d_name, "molecule_11181.sdf") == 0) {
            //printf("Ouverture du fichier %18s\n", entry->d_name);

            // Conversion en graphe moléculaire
            struct g_mol *molecule = gmol(file);
            
            // Sauvegarde du graphe moléculaire
            write_g_mol(molecule);
            //printf("Fichier %18s sauvegardé dans g_mol/g_mol_%d.dat \n", entry->d_name, molecule->Id);/*
            char filename[50];
            sprintf(filename, "g_mol/g_mol_%d.dat", molecule->Id);
            read_g_mol(filename);
            
            
            // Algorithme de McKay
            if (molecule->nb_atomes < 600) {
                SG_DECL(cg);
                mcKay(molecule, &cg);
                //print_sg(&cg);

                int aie = 0;
                for (size_t i = 0; i < cg.nde; i++) {
                    if (cg.e[i] >= cg.nv) {
                        printf("Mais c'est n'importe quoi !");
                        aie = 1;
                    }
                }

                // Plus courts chemins
                if (aie == 0) {
                    chemin **tab = malloc(molecule->nb_atomes * sizeof(chemin *));
                    for (int i = 0; i < molecule->nb_atomes; i++) {
                        tab[i] = malloc(molecule->nb_atomes * sizeof(chemin));
                        for (int j = 0; j < molecule->nb_atomes; j++) {
                            init_path(&tab[i][j], molecule->nb_atomes,
                                    2 * molecule->nb_liaisons);
                        }
                    }

                    smallest_paths(&cg, tab);

                    struct Cycle *Base;
                    int tailleB;

                    Base = Horton(&cg, tab, &tailleB);

                    //printf("\n\nFLAGGYFLAGGYFLAGFLAG FIN HORTON %d\n\n", tailleB);
                    //struct g_cycles * GrapheCycle;

                    TabGrapheCycle[nb_mol_vues] = ConvertBaseIntoGraph(molecule, Base, tailleB, &cg);

                    for (int i = molecule->nb_atomes - 1; i >= 0; i--) {
                        for (int j = molecule->nb_atomes - 1; j >= 0; j--) {
                            free_path(&tab[i][j]);
                        }
                        free(tab[i]);
                    }
                    free(tab);
                }
            
                SG_FREE(cg);
            }

            // Free et fermeture du fichier
            freeMolecule(molecule, molecule->nb_atomes);
            //printf("Fermeture du fichier %18s\n", entry->d_name);
            // }
            fclose(file);
        }

    }



/*
 ____________________________________________________________________________________
|                                                                                    |
|                COMPARAISON              ENTRE               MOLÉCULES              |
|____________________________________________________________________________________|

*/

calculSimilarite(TabGrapheCycle, nb_mol_vues);
graphesSimilaires(*TabGrapheCycle[5]);
/*
printf("\n\nFLAGGY FLAGGY FLAG FLAG COMPARAISON %d\n\n", nb_mol_vues);

    int id1 = 4;
    int id2 = 4;

    for (int i = 3; i < nb_mol_vues; i++) {
      // printf("\n\nVICTOIRE MOUSAILLON\n\n");

      if(TabGrapheCycle[i]->Id == 1296) {
        id1 = i;
      }

      if(TabGrapheCycle[i]->Id == 1712) {
        id2 = i;
      }
    }

    bool test = false;

    test = comparaisonNbreDeCycles(*TabGrapheCycle[id1], *TabGrapheCycle[id2]);


    if(test == true){
      printf("\n\nLES MOLÉCULES %d ET %d SONT IDENTIQUES !! VICTOIRE !!\n\n",TabGrapheCycle[id1]->Id,TabGrapheCycle[id2]->Id);
    } else {
      printf("\n\nBOUHOUHOUHÇAMARCHEPAS :/:/:/:/:/:/:/:/:/:/:/:/:/:/:/:/:/:/\n\n");

    }*/

    //classeEquivalences(TabGrapheCycle,nb_mol_vues);
    closedir(directory);
    return 0;
}