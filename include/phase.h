#ifndef PHASE_H
#define PHASE_H

#include "graph.h"
#include "syncmer.h"

#include <vector>

const int UNASSIGNED = -1;
const int DISCARDED = -2;
const int HET = 0;
const int HAP1 = 1;
const int HAP2 = 2;

namespace phase{
    void assign_initial_node_partitions(Graph* graph, std::vector<SketchedSeq*>& reads);
    bool reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads);
}

#endif
