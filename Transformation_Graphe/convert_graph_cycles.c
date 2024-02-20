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

      if(verifIntersectionVide(&c[i][A[j].IdA1], &c[i][A[j].IdA2], sg->nde)) {

        cpt++;

      }
    }
  }

  struct Cycle Ci[cpt];
  cpt = 0;

  for (int i = 0; i < sg->nv; i++) {
    for (int j = 0; j <nbA; j++) {

      if(i != A[j].IdA1 && i != A[j].IdA2) {
        if(verifIntersectionVide(&c[i][A[j].IdA1], &c[i][A[j].IdA2], sg->nde)) {
          
          // printf("sommet %d\t atome 1 : %d\t atome 2: %d\n",i,A[j].IdA1,A[j].IdA2);
          ajoutCycles(Ci, &c[i][A[j].IdA1], &c[i][A[j].IdA2], A[j], cpt, sg->nde, sg);
          cpt++;

        }
      }
    }
  }
  

  TriCroissant(Ci, cpt);

  /*printf("\n\nFLAGGY FLAG APRÈS TRI\n\n");

  for(int i = 0; i<cpt; i++) {
    printf("Cycle numéro %d, taille : %d\n",i,Ci[i].taille);
  }*/

  struct Cycle *Base;
  Base = malloc(cpt * sizeof(struct Cycle));
  if (Base == NULL) {
    fprintf(stderr, "Erreur : échec de l'allocation mémoire.\n");
    exit(EXIT_FAILURE);
}

  int tailleB = ExtractionBase(Ci, sg->nde, Base, cpt);

  /*printf("\n\nFLAGGY FLAG APRÈS EXTRACTION BASE\n\n");

  for(int i = 0; i<tailleB; i++) {
     printf("Cycle numéro %d, taille : %d\n",Base[i].Id,Base[i].taille);
         for (int j = 0; j < sg->nde; j++) {
        //printf("%d, ", Base[i].liaisons[j]);
    }
  }*/

  printf("\nTaille de la base de cycle : %d\t",tailleB);
  printf("\ttaille des cycles de la base : ");
  for(int i = 0; i < tailleB; i++){
    printf(" %d ,",Base[i].taille);
  }
  printf("\n");

}

  void init_liaisons(sparsegraph *sg, struct Liaison A[]) {

    int cptnbArete = 0;
    int cptIDA1 = 0;
    while(sg->d[cptIDA1] == 0) {
      cptIDA1++;
    }
    int *TabOccurence = malloc(sg->nv * sizeof(int));

    if (TabOccurence == NULL) {
      // Gérer l'échec de l'allocation mémoire
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sg->nv; i++) {
        TabOccurence[i] = 0;
    }


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

    free(TabOccurence);

    printf("\nNB LIAISONS = %d\n",cptnbArete);

    for(int i = 0; i < cptnbArete; i++) {
      //printf("liaison %d, atome 1 :%d\tatome 2 : %d\n",i, A[i].IdA1,A[i].IdA2);
    }
  }

int verifIntersectionVide(chemin *c1, chemin *c2, int m) {

  for (int i = 0; i < m; i++) {
    if ((c1->liaisons[i] == 1) && (c2->liaisons[i] == 1) ) {
      return 0;
    }
  }

  // printf("chemin :: atome 1 : %d\tatome 2 : %d\n",c1->atomesIds[0],c2->atomesIds[0]);

  return 1;
}

void ajoutCycles(struct Cycle *Ci, chemin *c1, chemin *c2,struct Liaison a, int compteur, int m, sparsegraph *sg) {

  Ci[compteur].liaisons = malloc(m * sizeof(int));
  Ci[compteur].Id = compteur;
  Ci[compteur].taille = c1->taille + c2->taille + 1;
  Ci[compteur].degre = 0;



  for (int i = 0; i < m; i++) {
    if (c1->liaisons[i] == 1 || c2->liaisons[i] == 1) {
      Ci[compteur].liaisons[i] = 1;
    }
  }

  int indice = 0;

  for (int i = 0; i < sg->nv; i++) {
    for(int j = 0; j < sg->d[i]; j++) {

      if (Ci[compteur].liaisons[indice] == 1) {

        int IDA1 = i;
        int IDA2 = sg->e[indice];

        int indice2 = 0;

        for (int iprime = 0; iprime < sg->nv; iprime++) {
          for(int jprime = 0; jprime < sg->d[iprime]; jprime++) {

            if ((iprime == IDA2) && (sg->e[indice2] == IDA1)) {

              Ci[compteur].liaisons[indice2] = 1;

            }
            indice2++;
          }
        }
      }


    if ((i == a.IdA1) && (sg->e[indice] == a.IdA2)) {
      Ci[compteur].liaisons[indice] = 1;

    }

    if ((i == a.IdA2) && (sg->e[indice] == a.IdA1)) {
      Ci[compteur].liaisons[indice] = 1;
    }

    indice++;
   }
  }

  //printf("\nLe cycle à bien été ajouté\n");

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

int  ExtractionBase(struct Cycle *Ci, int m, struct Cycle *Base, int tailleCi) {

  int *incidence = malloc(m * sizeof(int));
  int tailleBase = 0;

    if (incidence == NULL) {
      // Gérer l'échec de l'allocation mémoire
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < m; i++) {
        incidence[i] = 0;
    }


    for (int i = 0; i < tailleCi; i++) {

      int oui = 0;


      for (int j = 0; j < m; j++) {
        if (Ci[i].liaisons[j] == 1 && incidence[j] == 0) {
          oui = 1;
        }
      }

      if (oui == 1) {

        Base[tailleBase] = Ci[i];

        for (int j = 0; j < m; j++) {
          if (Base[tailleBase].liaisons[j] == 1 ) {
            incidence[j] = 1;
          }
        }

        tailleBase++;

      }

    }


    free(incidence);

    return tailleBase;
}