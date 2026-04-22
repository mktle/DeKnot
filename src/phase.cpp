#include "phase.h"

#include <cassert>
#include <iostream>

const int MAX_BACKTRACKS = 5;

void phase::initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs){
    graph->update_median_coverage();
    /************************
     * Map reads to unitigs *
     ************************/
    for (size_t i{0}; i < reads.size(); i++){
        if (reads[i]->is_read && reads[i]->kmers.size() >= (size_t)graph->get_k_prime()){
            ankerl::unordered_dense::set<Unitig*> alned_to;
            Unitig* track_utg {NULL};
            
            size_t check_thresh = reads[i]->kmers.size() - (size_t)graph->get_k_prime() + 1;
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};
                Node* curr_node = graph->get_node_from_lookup_table(curr_canon);
                if (curr_node != NULL && track_utg != curr_node->get_unitig()){
                    alned_to.insert(curr_node->get_unitig());
                    track_utg = curr_node->get_unitig();
                }
            }
            // now add read record to all unitigs
            for (const auto& aln_utg : alned_to) {
                aln_utg->add_read(reads[i]);
            }
        }
    }
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL){
            utg->set_assignment(UNASSIGNED);
            utg->set_parallel_utg();
            sorted_unitigs.push_back(utg);
        }
    }
    std::sort(sorted_unitigs.begin(), sorted_unitigs.end(), compare_unitigs);
}

bool phase::seed_initial_haplotypes(std::vector<Unitig*>& sorted_unitigs){
    /***************************************
     * Choose bubble to start phasing from *
     ***************************************/
    Unitig* initial_utg = NULL;
    Unitig* parallel_utg = NULL;
    for (size_t i{0}; i < sorted_unitigs.size(); i++){
        parallel_utg = sorted_unitigs[i]->get_parallel_utg();
        if (parallel_utg != NULL){
            initial_utg = sorted_unitigs[i];
            break;
        }
    }
    if (initial_utg == NULL) return false;
    initial_utg->set_assignment(HAP1);
    initial_utg->assign_reads(true);
    initial_utg->lock_haplotype();
    initial_utg->tally_reads();
    parallel_utg->set_assignment(HAP2);
    parallel_utg->assign_reads(true);
    parallel_utg->lock_haplotype();
    parallel_utg->tally_reads();
    return true;
}

// return true if any assignments changed
bool phase::reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads){
    /**
     * Clear the 
     *******/
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        if (unitig_list[i] != NULL){
            unitig_list[i]->set_hap_flipped(false);
        }
    }
     
    /**************************
     * Reassign unitig labels *
     **************************/
    bool changed {false};
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL){
            bool flipped = graph->reassign_haplotype(utg);
            if (flipped) changed = true;
        }
    }

    /************************
     * Reassign read labels *
     ************************/
    for (size_t i{0}; i < reads.size(); i++){
        if (reads[i]->is_read && !reads[i]->locked && reads[i]->kmers.size() >= (size_t)graph->get_k_prime()){
            int num1 {0};
            int num2 {0};
            size_t check_thresh = reads[i]->kmers.size() - (size_t)graph->get_k_prime() + 1;
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};
                Node* curr_node = graph->get_node_from_lookup_table(curr_canon);
                if (curr_node != NULL && curr_node->get_unitig()->get_parallel_utg() != NULL){
                    if (curr_node->get_unitig()->get_assignment() == HAP1) num1++;
                    if (curr_node->get_unitig()->get_assignment() == HAP2) num2++;
                }
            }
            if (num1 > num2 && (reads[i]->assignment == HAP2 || reads[i]->assignment == UNASSIGNED)){
                reads[i]->assignment = HAP1;
                changed =  true;
            }else if (num2 > num1 && (reads[i]->assignment == HAP1 || reads[i]->assignment == UNASSIGNED)){
                reads[i]->assignment = HAP2;
                changed = true;
            }
        }
    }
    return changed;
}

void phase::remove_tagged_unitigs(std::vector<Unitig*>& utg_list){
    auto assigned = [](Unitig* utg) { return utg->get_assignment() != UNASSIGNED; };
    auto new_end = std::remove_if(utg_list.begin(), utg_list.end(), assigned);
    utg_list.erase(new_end, utg_list.end());
}

void phase::get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs){
    /***************************************************
     * Identify which unitigs are at the end of a walk *
     ***************************************************/
    std::vector<Unitig*> hap1_endpoints;
    std::vector<Unitig*> hap2_endpoints;
    std::vector<Unitig*> hom_endpoints;

    const std::vector<Unitig*> utgs = graph->view_unitigs();
    for (size_t i{0}; i < utgs.size(); i++){
        if(utgs[i] != NULL){
            assert(utgs[i]->get_front() != utgs[i]->get_back());
            // subtract one for the neighbor on the same unitig
            int front_neighs {utgs[i]->get_front()->num_neighbors() - 1};
            int back_neighs {utgs[i]->get_back()->num_neighbors() - 1};

            if ((front_neighs == 0 && back_neighs > 0) || (back_neighs == 0 && front_neighs > 0)){
                if (utgs[i]->get_assignment() == HOM){
                    hom_endpoints.push_back(utgs[i]);
                }else if (utgs[i]->get_assignment() == HAP1){
                    hap1_endpoints.push_back(utgs[i]);
                }else if (utgs[i]->get_assignment() == HAP2){
                    hap2_endpoints.push_back(utgs[i]);
                }
            }
        }
    }
    std::sort(hap1_endpoints.begin(), hap1_endpoints.end(), compare_unitigs);
    std::sort(hap2_endpoints.begin(), hap2_endpoints.end(), compare_unitigs);
    std::sort(hom_endpoints.begin(), hom_endpoints.end(), compare_unitigs);

    std::vector<std::vector<Unitig*>*> hap_endpoints = {&hap1_endpoints, &hap2_endpoints};

    /**********************
     * Get walk sequences *
     **********************/
    graph->update_median_coverage();
     for (size_t i = 0; i < 2; i++) {
        graph->reset_haplotypes();
        auto& endpoints = *hap_endpoints[i];
        auto& other_endpoints = *hap_endpoints[1 - i];

        int hap_num = (i == 0) ? HAP1 : HAP2;
        int other_hap_num = (i == 0) ? HAP2 : HAP1;
        std::string walk_seq;
        Unitig* curr_utg = NULL;
        // if none in this hap, both ends must be homozygous
        if (endpoints.size() > 0){
            curr_utg = endpoints[0];
        }else if (hom_endpoints.size() > 0){
            curr_utg = hom_endpoints[0];
        }else{
            // navigate to first 
            Unitig* start_utg = other_endpoints[other_endpoints.size() - 1];
            Node* start_node = start_utg->get_front()->neighbors_diff_utg() ? start_utg->get_front() : start_utg->get_back();
            const std::vector<Edge>& neighs = start_node->view_neighbors();
            for (size_t j{0}; j < neighs.size(); j++){
                Unitig* neigh_utg = neighs[j].neighbor->get_unitig();
                if (neigh_utg != start_utg && (neigh_utg->get_assignment() == hap_num || (neigh_utg->get_assignment() == HOM && curr_utg == NULL))){
                    curr_utg = neigh_utg;
                }
            }
        }
        Node* opp_node = curr_utg->get_front()->neighbors_diff_utg() ? curr_utg->get_front() : curr_utg->get_back();
        curr_utg->get_hpc_sequence(walk_seq, opp_node == curr_utg->get_front(), false);
        // while there are still more unitigs along the walk
        while(opp_node->num_neighbors() > 1){
            Node* next_end {NULL};
            const std::vector<Edge>& opp_neighbors = opp_node->view_neighbors();
            std::vector<Node*> hap_neighs;
            std::vector<Node*> hom_neighs;
            std::vector<Node*> other_neighs;
            // find neighbor haplotype assignments
            for (size_t i{0}; i < opp_neighbors.size(); i++){
                Node* neigh = opp_neighbors[i].neighbor;
                Unitig* neigh_utg = neigh->get_unitig();
                if (neigh_utg->get_id() != opp_node->get_unitig()->get_id()){
                    if (neigh_utg->get_assignment() == hap_num) hap_neighs.push_back(neigh);
                    if(neigh_utg->get_assignment() == HOM) hom_neighs.push_back(neigh);
                    if(neigh_utg->get_assignment() == other_hap_num) other_neighs.push_back(neigh);
                }
            }
            // prioritize paths that are this haplotype, then homozygous, then other if needed
            std::vector<Node*>* choices = NULL;
            if (!hap_neighs.empty()) {
                choices = &hap_neighs;
            } else if (!hom_neighs.empty()) {
                choices = &hom_neighs;
            } else if (!other_neighs.empty()) {
                choices = &other_neighs;
            }
            assert(choices != NULL);
            if(choices->size() == 1){
                next_end = (*choices)[0];
            }else{
                int best_score = -1;
                Node* best_choice = NULL;
                float best_frac = 0;
                int best_other_num = 0;
                for (Node* n : *choices){
                    Unitig* n_utg = n->get_unitig();
                    Node* curr_op = n_utg->opposite_end(n);
                    int times_walked = n_utg->get_times_walked();
                    if (times_walked < n_utg->get_walk_limit()){
                        n_utg->set_times_walked(times_walked + 1);
                        int curr_score = get_total_score(curr_op, hap_num, 0, MAX_BACKTRACKS);
                        float curr_frac = n_utg->get_hap_frac(hap_num);
                        int other_num = n_utg->get_hap_count(other_hap_num);
                        if ((curr_score > best_score) || ((curr_score == best_score) && (curr_frac > best_frac)) || ((curr_score == best_score) && (curr_frac == best_frac) && (other_num < best_other_num))){
                            best_score = curr_score;
                            best_choice = n;
                            best_frac = curr_frac;
                            best_other_num = other_num;
                        }
                    }
                }
                next_end = best_choice;
            }
            if (next_end == NULL) break;
            next_end->get_unitig()->set_times_walked(next_end->get_unitig()->get_times_walked() + 1);
            next_end->get_unitig()->get_hpc_sequence(walk_seq, next_end == next_end->get_unitig()->get_back(), true);
            opp_node = next_end->get_unitig()->opposite_end(next_end);
        }
        ret_seqs.emplace_back(walk_seq);
    }
}

// called recursively to identify best local choce when multiple options
int phase::get_total_score(Node* start_node, int hap, int start_score, int steps_left){
    if (steps_left == 0){
        return start_score;
    }
    int best_score = start_score;
    std::vector<Edge> neighs = start_node->view_neighbors();
    for (size_t i{0}; i < neighs.size(); i++){
        Node* neigh_node = neighs[i].neighbor;
        Unitig* neigh_utg = neigh_node->get_unitig();
        int times_walked = neigh_utg->get_times_walked();
        if (neigh_utg != start_node->get_unitig() && times_walked + 1 < neigh_utg->get_walk_limit()){
            int curr_weight = (neigh_utg->get_assignment() == hap || neigh_utg->get_assignment() == HOM) ? neigh_utg->num_nodes() : 0;
            neigh_utg->set_times_walked(times_walked + 1);
            Node* opp_node = neigh_utg->opposite_end(neigh_node);
            neigh_utg->set_times_walked(times_walked);
            int next_score = curr_weight + get_total_score(opp_node, hap, start_score + curr_weight, steps_left - 1);
            if (next_score > best_score){
                best_score = next_score;
            }
        }
    }
    return best_score;
}
