#include "convert_graph_cycles.h"
#include "gmol.h"
#include "comparaison.h"
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

void affiche_time(clock_t current_time) {

  clock_t test = clock();
  test = test - current_time;
  printf("%ld min %ld sec %ld ms\n",test/(CLOCKS_PER_SEC*60), (test/CLOCKS_PER_SEC)%60, (test*1000/CLOCKS_PER_SEC)%1000 );
  //printf("Finished in %ld µs\n\n",(test*1000000/CLOCKS_PER_SEC) );
}


int main(int argc, char *argv[]) {

  clock_t clock_start = clock(); 

    FILE *parametre = fopen(argv[1], "r");

    int mode = -1;
    int error = fscanf(parametre, "%d", &mode);
    if (error == EOF){
      printf("\nerreur scan\n");
      exit(1);
    }
    int IDmol1 = -1;
    error = fscanf(parametre, "%d", &IDmol1);
    if (error == EOF){
      printf("\nerreur scan\n");
      exit(1);
    }
    int IDmol2 = -1;
    error = fscanf(parametre, "%d", &IDmol2);
    if (error == EOF){
      printf("\nerreur scan\n");
      exit(1);
    }
    
    if (mode != 1 && mode != 2 && mode != 3) {
      printf("\n\nLe mode d'éxecution %d n'éxiste pas !!\n\n",mode);
      exit(10);
    }

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

  if(mode == 1) {
    TabGrapheCycle = malloc (cpt * sizeof(struct g_cycles));
  } else if (mode == 2) {
    TabGrapheCycle = malloc (sizeof(struct g_cycles));
  } else if (mode == 3) {
    TabGrapheCycle = malloc (2 * sizeof(struct g_cycles));
  }



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

    while ((entry = readdir(directory)) != NULL) {
       /* printf("--------------------------------- Molécule n°%d "
                "---------------------------------\n",
                nb_mol_vues);*/
        nb_mol_vues++;

        // Vérifie si le fichier est un fichier .sdf
        if (strstr(entry->d_name, ".sdf") != NULL) {


         // if(mode == 1 || (mode == 2 && ) {
            // Construit le chemin complet du fichier
            char filepath[512];
            snprintf(filepath, sizeof(filepath),
                    "../TraitementDonnee/output_molecules_v3000/%s",
                    entry->d_name); // à modifier quand ce sera bon
            FILE *file = fopen(filepath, "r");
            if (!file) {
                printf("Impossible d'ouvrir le fichier");
            }

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

              if((mode == 1) || (mode == 2 && molecule->Id == IDmol1) || (mode == 3 && molecule->Id == IDmol1) || (mode == 3 && molecule->Id == IDmol2)) {
                
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

                      int indiceTab = nb_mol_vues;
                      if(mode == 2) {
                        indiceTab = 0;
                      }
                      if(mode == 3 && molecule->Id == IDmol1) {
                        indiceTab = 0;

                      }
                      if(mode == 3 && molecule->Id == IDmol2) {
                        indiceTab = 1;

                      }

                      TabGrapheCycle[indiceTab] = ConvertBaseIntoGraph(molecule, Base, tailleB, &cg);

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
            }

            // Free et fermeture du fichier
            freeMolecule(molecule, molecule->nb_atomes);
            //printf("Fermeture du fichier %18s\n", entry->d_name);
            fclose(file);
            
          
        }

    }



/*
 ____________________________________________________________________________________
|                                                                                    |
|                COMPARAISON              ENTRE               MOLÉCULES              |
|____________________________________________________________________________________|

*/

    if (mode == 1) {
    printf("\n\nCONVERSION EN GRAPHE DE CYCLE TERMINÉ - DÉBUT CALCUL SIMILARITÉ\n");

    printf("durée : ");
    affiche_time(clock_start);

    calculSimilarite(TabGrapheCycle, nb_mol_vues);

    printf("\n\nCALCUL SIMILARITÉ TERMINÉ - DÉBUT CALCUL CLASSES D'ÉQUIVALENCES\n");

    printf("durée : ");
    affiche_time(clock_start);

    classeEquivalences(TabGrapheCycle,nb_mol_vues);
    printf("\n\nCALCUL CLASSES D'ÉQUIVALENCES TERMINÉ\n\n");
    }

    if(mode == 2) {
      graphesSimilaires(*TabGrapheCycle[0]);
    }

    if(mode == 3) {
      double score = 0.0;
      score = similarite(*TabGrapheCycle[0],*TabGrapheCycle[1]);
      printf("\n\nLes molécules d'ID chEBI %d et %d ont un score de similarité de : %lf\n\n",TabGrapheCycle[0]->Id,TabGrapheCycle[1]->Id,score);
    }

    closedir(directory);

    printf("\n\nFin du programme\n\n");
    printf("durée : ");
    affiche_time(clock_start);
    return 0;
}