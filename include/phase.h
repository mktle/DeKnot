#ifndef PHASE_H
#define PHASE_H

#include "graph.h"
#include "syncmer.h"

#include <vector>

const int UNASSIGNED = -1;
const int HOM = 0;
const int HAP1 = 1;
const int HAP2 = 2;

struct Walk{
    Walk* prefix;
    std::vector<Unitig*> utgs;
    int num_hap_nodes;

    Walk(Walk* p, int n): prefix(p), utgs(), num_hap_nodes(n) {}
};

namespace phase{
    void initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs);
    bool seed_initial_haplotypes(std::vector<Unitig*>& sorted_unitigs);
    bool reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads);
    void remove_tagged_unitigs(std::vector<Unitig*>& utg_list);
    void identify_endpoints(Graph* graph, std::vector<std::vector<Unitig*>>& hap_endpoints, std::vector<Unitig*>& hom_endpoints);
    void get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs);
    void add_to_walk(std::vector<Walk>& all_walks, Walk* curr_walk, Node* to_add, int hap_num);
}

#endif
