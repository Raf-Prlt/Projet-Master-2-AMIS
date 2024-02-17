#include "convert_graph_cycles.h"
#include "gmol.h"
#include "nauty2_8_8/nausparse.h" /* which includes nauty.h */

DYNALLSTAT(int, lab, lab_sz);
DYNALLSTAT(int, ptn, ptn_sz);
DYNALLSTAT(int, orbits, orbits_sz);
static DEFAULTOPTIONS_SPARSEGRAPH(options);

void gmol_to_sparsegraph(struct g_mol *g, sparsegraph *sg) {
    sg->nv = g->nb_atomes;
    sg->nde = 2 * g->nb_liaisons;
    int j = 0;
    for (int i = 0; i < g->nb_atomes; i++) {
        struct Atome a = g->Atomes[i];
        sg->d[i] = a.degre;
        sg->v[i] = j;
        for (int k = 0; k < a.degre; k++) {
            sg->e[j + k] = a.voisinsIds[k] - 1;
            // ici on doit avoir la propriété Atome[i].Id == i
        }
        j = j + a.degre;
    }
    print_sg(sg);
}

void print_sg(sparsegraph *sg) {
    printf("sg->nv = %d \n", sg->nv);
    printf("sg->nde = %ld \n", sg->nde);
    printf("sg->d = [");
    for (int i = 0; i < sg->nv; i++) {
        printf("%d, ", sg->d[i]);
    }
    printf("]\n");
    printf("sg->v = [");
    for (int i = 0; i < sg->nv; i++) {
        printf("%ld, ", sg->v[i]);
    }
    printf("]\n");
    printf("sg->e = [");
    for (int i = 0; i < sg->nde; i++) {
        printf("%d, ", sg->e[i]);
    }
    printf("]\n");
}

void print_path(chemin *c, int m) {
    printf("taille = %d, ", c->taille);
    printf("atomesIds = [");
    for (int i = 0; i < c->taille; i++) {
        printf("%d, ", c->atomesIds[i]);
    }
    printf("]\n");
    printf("liaisons = [");
    for (int i = 0; i < m; i++) {
        printf("%d, ", c->liaisons[i]);
    }
    printf("]\n");
}

void mcKay(struct g_mol *g, sparsegraph *cg) {
    statsblk stats;
    sparsegraph sg;

    SG_INIT(sg);
    SG_ALLOC(sg, g->nb_atomes, 2 * g->nb_liaisons, "malloc");

    gmol_to_sparsegraph(g, &sg);

    options.getcanon = TRUE;

    DYNALLOC1(int, lab, lab_sz, g->nb_atomes, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, sg.nv, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, sg.nv, "malloc");

    sparsenauty(&sg, lab, ptn, orbits, &options, &stats, cg);

    SG_FREE(sg);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    DYNFREE(orbits, orbits_sz);
}

void init_path(chemin *c, int n, int m) {
    c->taille = 0;
    c->atomesIds = malloc(n * sizeof(int));
    c->liaisons = malloc(m * sizeof(int));
    for (int i = 0; i < m; i++) {
        c->liaisons[i] = 0;
    }
}

void free_path(chemin *c) {
    free(c->atomesIds);
    free(c->liaisons);
}

void smallest_paths(sparsegraph *sg, chemin **ps_cts_chms) {
    int n = sg->nv;
    int m = sg->nde;

    int *deja_vu = malloc(n * sizeof(int));

    int size_next = 0;
    int *next = malloc(n * sizeof(int));
    int size_current = 0;
    int *current = malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        // printf("-------i: %d---------\n", i);
        for (int j = 0; j < n; j++) {
            deja_vu[j] = 0;
        };
        deja_vu[i] = 1;
        next[0] = i;
        size_next = 1;

        ps_cts_chms[i][i].taille = 0;

        while (size_next > 0) {
            size_current = size_next;
            for (int j = 0; j < size_next; j++) {
                current[j] = next[j];
            }
            size_next = 0;
            for (int j = 0; j < size_current; j++) {
                int s = current[j];
                deja_vu[s] = 1;
                int degre = sg->d[s];
                int start = sg->v[s];
                int taille_chemin = ps_cts_chms[i][s].taille;
                int *atomesIds_chemin = ps_cts_chms[i][s].atomesIds;
                int *liaisons_chemin = ps_cts_chms[i][s].liaisons;
                for (int k = 0; k < degre; k++) {
                    int succ = sg->e[start + k];
                    // printf("succ: %d ", succ);
                    if (deja_vu[succ] == 0) {
                        // on ajoute succ à next s'il n'est pas déjà présent
                        int in = 0;
                        for (int l = 0; l < size_next; l++) {
                            if (next[l] == succ) {
                                in = 1;
                            }
                        }
                        if (in == 0) {
                            next[size_next] = succ;
                            size_next++;
                        }
                        chemin *c = &ps_cts_chms[i][succ];
                        if (c->taille == 0) {
                            c->taille = taille_chemin + 1;
                            for (int l = 0; l < taille_chemin; l++) {
                                c->atomesIds[l] = atomesIds_chemin[l];
                            }
                            for (int l = 0; l < m; l++) {
                                c->liaisons[l] = liaisons_chemin[l];
                            }
                            c->atomesIds[taille_chemin] = s;
                            c->liaisons[start + k] = 1;
                            // printf("(%d, %d) ", i, succ);
                            // print_path(c, m);
                            // printf("\n");
                            // fflush(stdout);
                        }
                    }
                }
            }
        }
    }

    free(deja_vu);
    free(current);
    free(next);
}

void Horton(sparsegraph *sg, chemin **c) {

  printf("\nHello there\nGENERAL KENOBI\n");

  printf("\nNb atome : %d\n",sg->nv);
  printf("\nNb liaisons : %ld\n",sg->nde);

  int nbA = sg->nde/2;

  struct Liaison A[nbA];

  init_liaisons(sg, A);

  int cpt= 0;
  

  for (int i = 0; i < sg->nv; i++) {
    for (int j = 0; j <nbA; j++) {
      printf("\nFLAGGYFLAGGYFLAGFLAG %d\n", cpt);

      if(verifIntersectionVide(&c[i][A[j].IdA1], &c[i][A[j].IdA2])) {

        cpt++;

      }
    }
  }

  struct Cycle Ci[cpt];
  cpt = 0;

  for (int i = 0; i < sg->nv; i++) {
    for (int j = 0; j <sg->nde/2; j++) {

      if(verifIntersectionVide(&c[i][A[j].IdA1], &c[i][A[j].IdA2])) {

        ajoutCycles(Ci, &c[i][A[j].IdA1], &c[i][A[j].IdA2], A[j], cpt);
        cpt++;

      }
    }
  }

  printf("\n\nFLAGGY FLAG AVANT TRI\n\n");

    for(int i = 0; i<cpt; i++) {
    printf("Cycle numéro %d, taille : %d\n",i,Ci[i].taille);
  }



  TriCroissant(Ci, cpt);

  printf("\n\nFLAGGY FLAG APRÈS TRI\n\n");


  for(int i = 0; i<cpt; i++) {
    printf("Cycle numéro %d, taille : %d\n",i,Ci[i].taille);
  }



}

  void init_liaisons(sparsegraph *sg, struct Liaison A[]) {

    int cptnbArete = 0;
    int cptIDA1 = 0;
    int TabOccurence[32] = {0};

    for(int i = 0; i < sg->nde; i++) {

      if (i >= sg->v[cptIDA1+1]) {
        cptIDA1++;
      }

      if (TabOccurence[cptIDA1] < sg->d[cptIDA1] && TabOccurence[sg->e[i]] < sg->d[sg->e[i]]) {

        A[cptnbArete].IdA1 = cptIDA1;
        A[cptnbArete].IdA2 = sg->e[i];
        A[cptnbArete].Poids = 1;

        cptnbArete++;
        TabOccurence[cptIDA1]++;
        TabOccurence[sg->e[i]]++;
      }

    }

    printf("\nNB LIAISONS = %d\n",cptnbArete);
  }

int verifIntersectionVide(chemin *c1, chemin *c2) {

  for (int i = 1; i < c1->taille; i++) {
    for (int j = 1; j < c2->taille; j++) {

        if(c1->atomesIds[i] == c2->atomesIds[j]) {
          return 0;
        }
    }
  }

  return 1;
}

void ajoutCycles(struct Cycle *Ci, chemin *c1, chemin *c2,struct Liaison a, int compteur) {

  Ci[compteur].liaisons = malloc (100 * sizeof(int));

  Ci[compteur].Id = compteur;
  Ci[compteur].taille = c1->taille + c2->taille + 1;
  Ci[compteur].degre = 0;

  int cpt = 0;

  for (int i = 0; i < c1->taille; i++) {
    Ci[compteur].liaisons[cpt] = c1->liaisons[i];
    cpt++;
  }

  for (int i = 0; i < c2->taille; i++) {
    Ci[compteur].liaisons[cpt] = c2->liaisons[i];
    cpt++;
  }

    //Ci[compteur].liaisons[cpt] = a;

  //printf("\nLe cycle à bien été ajouté\n");
  //printf("\nCompteur : %d\tverif : %d\n",compteur, 29*32);

}

void TriCroissant(struct Cycle Ci[], int tailleTab) {



  for(int i = 0; i < tailleTab; i++) {

    int tailleMin = 100000000;
    int indiceMin = 0;

    for(int j = i; j < tailleTab; j++) {

      if(Ci[j].taille < tailleMin) {
        indiceMin = j;
        tailleMin = Ci[j].taille;
      }

    }

    struct Cycle tmp;
    tmp = Ci[indiceMin];
    Ci[indiceMin] = Ci[i];
    Ci[i] = tmp;

  }

}