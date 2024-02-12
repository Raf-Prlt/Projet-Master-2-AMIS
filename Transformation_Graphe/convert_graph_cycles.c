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
      sg->e[j + k] =
          a.voisinsIds[k]; // ici on doit avoir la propriété Atome[i].Id == i
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
      // printf("size_current = %d \n", size_current);
      // printf("current = [");
      // for (int j = 0; j < size_current; j++) {
      //   printf("%d, ", current[j]);
      // }
      // printf("]\n");
      for (int j = 0; j < size_current; j++) {
        int s = current[j];
        // printf("s: %d, ", s);
        // fflush(stdout);
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
              // printf("oxo succ: %d, size_next: %d ", succ, size_next);
              next[size_next] = succ;
              size_next++;
            }
            chemin *c = &ps_cts_chms[i][succ];
            // printf("c->taille = %d, ", c->taille);
            if (c->taille == 0) {
              c->taille = taille_chemin + 1;
              for (int l = 0; l < taille_chemin; l++) {
                c->atomesIds[l] = atomesIds_chemin[l];
              }
              // printf("coucou7 ");
              // fflush(stdout);
              for (int l = 0; l < m; l++) {
                c->liaisons[l] = liaisons_chemin[l];
              }
              // printf("coucou8 ");
              // fflush(stdout);
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