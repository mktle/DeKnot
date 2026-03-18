#ifndef PHASE_H
#define PHASE_H

#include "graph.h"
#include "syncmer.h"

#include <vector>

const int UNASSIGNED = -1;
const int HOM = 0;
const int HAP1 = 1;
const int HAP2 = 2;

namespace phase{
    void initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs);
    void seed_initial_haplotypes(std::vector<Unitig*>& sorted_unitigs);
    bool reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads);
    void remove_tagged_unitigs(std::vector<Unitig*>& utg_list);
    void get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs);
}

#endif
