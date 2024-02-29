#include "convert_graph_cycles.h"
#include "gmol.h"
#include <stdbool.h>

double comparaisonNbreDeCycles(struct g_cycles g1, struct g_cycles g2);
double comparaisonTailleDeCycles(struct g_cycles g1, struct g_cycles g2);
double comparaisonNbreVoisinsDeCycle(struct g_cycles g1, struct g_cycles g2);
double comparaisonPoidsDesArretes(struct g_cycles g1, struct g_cycles g2);
double similarite(struct g_cycles g1, struct g_cycles g2);
void calculSimilarite(struct g_cycles **liste, int size);
void graphesSimilaires(struct g_cycles g1);
int rechercheIndiceAvecId(struct g_cycles **TabGrapheCycle, int cpt, int id);
void classeEquivalences(struct g_cycles **TabGrapheCycle, int cpt);