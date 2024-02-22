#include "convert_graph_cycles.h"
#include "gmol.h"
#include <stdbool.h>

bool comparaisonNbreDeCycles(struct g_cycles g1, struct g_cycles g2);
double comparaisonTailleDeCycles(struct g_cycles g1, struct g_cycles g2);
double comparaisonNbreVoisinsDeCycle(struct g_cycles g1, struct g_cycles g2);
double degreDeSimilarite(struct g_cycles g1, struct g_cycles g2);
void calculSimilarites(struct g_cycles *liste);
void similarites(struct g_cycles *liste, struct g_cycles g1);
int mainComp();