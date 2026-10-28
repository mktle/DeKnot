#ifndef PHASE_H
#define PHASE_H

#include "graph.h"
#include "syncmer.h"

#include <vector>

const int LOOKAHEAD_BUDGET = 256;

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
    void initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs, int comp_num);
    bool seed_initial_haplotypes(std::vector<Unitig*>& sorted_unitigs);
    bool reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& utg_list);
    void remove_tagged_unitigs(std::vector<Unitig*>& utg_list);
    void identify_endpoints(Graph* graph, std::vector<std::vector<Unitig*>>& hap_endpoints, std::vector<Unitig*>& hom_endpoints, int comp_num);
    bool get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs, int comp_num, std::vector<Unitig*>* walked = NULL);
    void build_walk(Walk& walk, Node* start_node, int hap_num);
    int lookahead_score(Node* to_add, int hap_num, int& budget);
    void phase_component(int comp_num, std::vector<std::string>& phased_walks, Graph* prev_graph, std::vector<SketchedSeq*>& corrected_seqs, bool verbose);
    void assign_reads_by_syncmers(Graph* graph, std::vector<SketchedSeq*>& reads, const std::vector<uint64_t>& neutral, bool verbose);
}

#endif
