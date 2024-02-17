#include "convert_graph_cycles.h"
#include "gmol.h"
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    DIR *directory;
    struct dirent *entry;
    directory = opendir(
        "../TraitementDonnee/output_molecules/"); // à modifier quand ce sera
                                                // bon (remplacer test par
                                                // output_molecules)
    if (!directory) {
    printf("Impossible d'ouvrir le dossier");
    exit(1);
    }
    int nb_mol_vues = 0;
    // Lit chaque fichier dans le dossier
    while ((entry = readdir(directory)) != NULL) {
    /* printf("--------------------------------- Molécule n°%d "
             "---------------------------------\n",
             nb_mol_vues);*/
    nb_mol_vues++;

    // Vérifie si le fichier est un fichier .sdf
    if (strstr(entry->d_name, ".sdf") != NULL) {

        // Construit le chemin complet du fichier
        char filepath[512];
        snprintf(filepath, sizeof(filepath),
                 "../TraitementDonnee/output_molecules/%s",
                 entry->d_name); // à modifier quand ce sera bon
        FILE *file = fopen(filepath, "r");
        if (!file) {
        printf("Impossible d'ouvrir le fichier");
        }

       if(strcmp(entry->d_name, "molecule_7.sdf") == 0) {


        printf("Ouverture du fichier %18s\n", entry->d_name);

        // Conversion en graphe moléculaire
        struct g_mol *molecule = gmol(file);

        // Sauvegarde du graphe moléculaire
        if ((molecule->Id != -1) && (molecule->name) && (molecule->nb_atomes) &&
            (molecule->nb_liaisons)) {
        write_g_mol(molecule);
        printf("Fichier %18s sauvegardé dans g_mol/g_mol_%d.dat \n",
                 entry->d_name, molecule->Id);
        char filename[50];
        sprintf(filename, "g_mol/g_mol_%d.dat", molecule->Id);
        read_g_mol(filename);
        }

        // Algorithme de McKay
        if (molecule->nb_atomes < 600) {
        SG_DECL(cg);
        mcKay(molecule, &cg);
        print_sg(&cg);

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

            Horton(&cg, tab);


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
        printf("Fermeture du fichier %18s\n", entry->d_name);
        freeMolecule(molecule, molecule->nb_atomes);
    }
        fclose(file);
    }
    }
    closedir(directory);
    return 0;
}