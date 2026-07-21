#include "phase.h"

#include <cassert>
#include <iostream>

void phase::initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs){
    graph->update_median_coverage();
    if (!graph->aligned_reads()) graph->align_reads_to_unitigs(reads);
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
                //if (curr_node != NULL && curr_node->get_unitig()->get_parallel_utg() != NULL){
                if (curr_node != NULL){
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

void phase::identify_endpoints(Graph* graph, std::vector<std::vector<Unitig*>>& hap_endpoints, std::vector<Unitig*>& hom_endpoints){
    hap_endpoints.resize(2);
    auto& hap1_endpoints = hap_endpoints[0];
    auto& hap2_endpoints = hap_endpoints[1];

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
}

void phase::get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs){
    std::vector<std::vector<Unitig*>> hap_endpoints;
    std::vector<Unitig*> hom_endpoints;
    phase::identify_endpoints(graph, hap_endpoints, hom_endpoints);

    /**********************
     * Get walk sequences *
     **********************/
    Node* start_node = NULL;
    graph->update_median_coverage();
     for (size_t i = 0; i < 2; i++) {
        graph->reset_haplotypes();
        auto& endpoints = hap_endpoints[i];
        auto& other_endpoints = hap_endpoints[1 - i];

        int hap_num = (i == 0) ? HAP1 : HAP2;
        // if none in this hap, both ends must be homozygous
        if (endpoints.size() > 0){
            start_node = endpoints[0]->get_front()->neighbors_diff_utg() ? endpoints[0]->get_back() : endpoints[0]->get_front();
        }else if (hom_endpoints.size() > 0){
            start_node = hom_endpoints[0]->get_front()->neighbors_diff_utg() ? hom_endpoints[0]->get_back() : hom_endpoints[0]->get_front();
        }else{
            // navigate to first
            Unitig* start_utg = other_endpoints[other_endpoints.size() - 1];
            Node* first_node = start_utg->get_front()->neighbors_diff_utg() ? start_utg->get_front() : start_utg->get_back();
            const std::vector<Edge>& neighs = first_node->view_neighbors();
            start_node = NULL;
            for (size_t j{0}; j < neighs.size(); j++){
                Node* neigh = neighs[j].neighbor;
                Unitig* neigh_utg = neigh->get_unitig();
                if (neigh_utg != start_utg && (neigh_utg->get_assignment() == hap_num || (neigh_utg->get_assignment() == HOM && start_node == NULL))){
                    start_node = neigh;
                    if (neigh_utg->get_assignment() == hap_num) break;
                }
            }
        }
        std::vector<Walk> all_walks;
        Walk* curr_walk = new Walk(NULL, 0);
        // recursive function that adds all walks
        add_to_walk(all_walks, curr_walk, start_node, hap_num);
        delete curr_walk;

        /*******************************************************
         * Identify best walk and concatenate unitig sequences *
         *******************************************************/
        Walk best_walk = all_walks[0];
        for (size_t j{1}; j < all_walks.size(); j++){
            if(all_walks[j].num_hap_nodes > best_walk.num_hap_nodes){
                best_walk = all_walks[j];
            }
        }
        std::string walk_seq;
        assert(best_walk.utgs[0] == start_node->get_unitig());
        start_node->get_unitig()->get_hpc_sequence(walk_seq, start_node == start_node->get_unitig()->get_back(), false);

        Node* exit_node = start_node->get_unitig()->opposite_end(start_node);
        for (size_t j{1}; j < best_walk.utgs.size(); j++){
            Node* next_node = (best_walk.utgs[j]->get_front()->is_neighbors_with(exit_node)) ? best_walk.utgs[j]->get_front() : best_walk.utgs[j]->get_back();
            next_node->get_unitig()->get_hpc_sequence(walk_seq, next_node == next_node->get_unitig()->get_back(), true);
            exit_node = next_node->get_unitig()->opposite_end(next_node);
        }
        ret_seqs.emplace_back(walk_seq);
    }
}

void phase::add_to_walk(std::vector<Walk>& all_walks, Walk* curr_walk, Node* to_add, int hap_num){
    int other_hap_num = (hap_num == HAP1) ? HAP2 : HAP1;
    Unitig* add_utg = to_add->get_unitig();
    int times_walked = add_utg->get_times_walked();
    // dead end
    if (times_walked >= add_utg->get_walk_limit()){
        all_walks.push_back(*curr_walk);
        return;
    }
    add_utg->set_times_walked(times_walked + 1);
    curr_walk->utgs.push_back(add_utg);
    int hap_weight = (add_utg->get_assignment() == hap_num || add_utg->get_assignment() == HOM) ? add_utg->num_nodes() : 0;
    curr_walk->num_hap_nodes += hap_weight;

    Node* opp_node = add_utg->opposite_end(to_add);
    const std::vector<Edge>& opp_neighbors = opp_node->view_neighbors();
    std::vector<Node*> hap_neighs;
    std::vector<Node*> other_neighs;
    // find neighbor haplotype assignments
    for (size_t i{0}; i < opp_neighbors.size(); i++){
        Node* neigh = opp_neighbors[i].neighbor;
        Unitig* neigh_utg = neigh->get_unitig();
        if (neigh_utg->get_id() != opp_node->get_unitig()->get_id()){
            if (neigh_utg->get_assignment() == hap_num) hap_neighs.push_back(neigh);
            if(neigh_utg->get_assignment() == HOM) hap_neighs.push_back(neigh);
            if(neigh_utg->get_assignment() == other_hap_num) other_neighs.push_back(neigh);
        }
    }
    // prioritize paths that are this haplotype or homozygous, then other if needed
    std::vector<Node*>* choices = NULL;
    if (!hap_neighs.empty()) {
        choices = &hap_neighs;
    } else if (!other_neighs.empty()) {
        choices = &other_neighs;
    }else{
        // no more neighbors to add
        all_walks.push_back(*curr_walk);
        add_utg->set_times_walked(times_walked);
        curr_walk->utgs.pop_back();
        curr_walk->num_hap_nodes -= hap_weight;
        return;
    }
    for (size_t i{0}; i < choices->size(); i++){
        Node* next = (*choices)[i];
        phase::add_to_walk(all_walks, curr_walk, next, hap_num);
    }
    add_utg->set_times_walked(times_walked);
    curr_walk->utgs.pop_back();
    curr_walk->num_hap_nodes -= hap_weight;
}
