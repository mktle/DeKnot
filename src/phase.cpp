#include "phase.h"

#include <cassert>
#include <iostream>
#include <unordered_set>

void phase::initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs, int comp_num){
    graph->update_median_coverage();
    if (!graph->aligned_reads()) graph->align_reads_to_unitigs(reads);
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL && utg->get_component() == comp_num){
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
bool phase::reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& utg_list){
    std::vector<int> before(utg_list.size());
    for (size_t i{0}; i < utg_list.size(); i++){
        utg_list[i]->set_hap_updated(false);
        before[i] = utg_list[i]->get_assignment();
    }

    /**************************
     * Reassign unitig labels *
     **************************/
    bool changed {false};
    for (size_t i{0}; i < utg_list.size(); i++){
        graph->reassign_haplotype(utg_list[i]);
    }
    for (size_t i{0}; i < utg_list.size(); i++){
        if (utg_list[i]->get_assignment() != before[i]){
            changed = true;
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
            utils::RollingHasher rh(reads[i]->kmers.data(), reads[i]->kmers.size(), (size_t)graph->get_k_prime());
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION, rh.hash_at(j, curr_orient == REVERSE_ORIENTATION)};
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

void phase::identify_endpoints(Graph* graph, std::vector<std::vector<Unitig*>>& hap_endpoints, std::vector<Unitig*>& hom_endpoints, int comp_num){
    hap_endpoints.resize(2);
    auto& hap1_endpoints = hap_endpoints[0];
    auto& hap2_endpoints = hap_endpoints[1];

    const std::vector<Unitig*> utgs = graph->view_unitigs();
    for (size_t i{0}; i < utgs.size(); i++){
        if(utgs[i] != NULL && utgs[i]->get_component() == comp_num){
            int front_neighs, back_neighs;
            if (utgs[i]->get_front() == utgs[i]->get_back()){
                int n_neighs {utgs[i]->get_front()->num_neighbors()};
                front_neighs = (n_neighs == 1) ? 0 : n_neighs;
                back_neighs = n_neighs;
            }else{
                // subtract one for the neighbor on the same unitig
                front_neighs = utgs[i]->get_front()->num_neighbors() - 1;
                back_neighs = utgs[i]->get_back()->num_neighbors() - 1;
            }

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

bool phase::get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs, int comp_num, std::vector<Unitig*>* walked){
    std::vector<std::vector<Unitig*>> hap_endpoints;
    std::vector<Unitig*> hom_endpoints;
    phase::identify_endpoints(graph, hap_endpoints, hom_endpoints, comp_num);
    if(hap_endpoints[0].size() == 0 && hap_endpoints[1].size() == 0 && hom_endpoints.size() == 0){
        return false;
    }

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
            const EdgeList& neighs = first_node->view_neighbors();
            start_node = NULL;
            for (size_t j{0}; j < neighs.size(); j++){
                Node* neigh = neighs[j].neighbor;
                Unitig* neigh_utg = neigh->get_unitig();
                if (neigh_utg != start_utg && (neigh_utg->get_assignment() == hap_num || (neigh_utg->get_assignment() == HOM && start_node == NULL))){
                    start_node = neigh;
                    if (neigh_utg->get_assignment() == hap_num) break;
                }
            }
            if (start_node == NULL){
                start_node = first_node;
            }
        }
        {
            std::vector<Node*> cands;
            for (Unitig* u : endpoints) cands.push_back(u->get_front()->neighbors_diff_utg() ? u->get_back() : u->get_front());
            for (Unitig* u : hom_endpoints) cands.push_back(u->get_front()->neighbors_diff_utg() ? u->get_back() : u->get_front());
            cands.push_back(start_node);
            long best_score {-1};
            Node* best_start = start_node;
            for (Node* cand : cands){
                graph->reset_haplotypes();
                Walk trial(NULL, 0);
                phase::build_walk(trial, cand, hap_num);
                long hap_nodes {0}, total_nodes {0};
                for (Unitig* u : trial.utgs){
                    total_nodes += u->num_nodes();
                    if (u->get_assignment() == hap_num) hap_nodes += u->num_nodes();
                }
                long score = hap_nodes * 1000000L + total_nodes;
                if (score > best_score){ best_score = score; best_start = cand; }
            }
            graph->reset_haplotypes();
            start_node = best_start;
        }
        Walk best_walk(NULL, 0);
        phase::build_walk(best_walk, start_node, hap_num);
        if (walked != NULL) walked->insert(walked->end(), best_walk.utgs.begin(), best_walk.utgs.end());

        std::string walk_seq;
        assert(best_walk.utgs[0] == start_node->get_unitig());
        auto single_node_reverse = [](Node* node, Node* other, bool entering) -> bool {
            const EdgeList& edges = entering ? other->view_neighbors() : node->view_neighbors();
            Node* target = entering ? node : other;
            for (size_t m{0}; m < edges.size(); m++){
                if (edges[m].neighbor == target){
                    return (entering ? edges[m].neighbor_orientation : edges[m].self_orientation) == REVERSE_ORIENTATION;
                }
            }
            return false;
        };
        Unitig* first_utg = start_node->get_unitig();
        bool first_reverse = (start_node == first_utg->get_back());
        if (first_utg->get_front() == first_utg->get_back() && best_walk.utgs.size() > 1){
            Unitig* second = best_walk.utgs[1];
            Node* second_entry = second->get_front()->is_neighbors_with(start_node) ? second->get_front() : second->get_back();
            first_reverse = single_node_reverse(start_node, second_entry, false);
        }
        first_utg->get_hpc_sequence(walk_seq, first_reverse, false);

        Node* exit_node = first_utg->opposite_end(start_node);
        for (size_t j{1}; j < best_walk.utgs.size(); j++){
            Unitig* next_utg = best_walk.utgs[j];
            Node* next_node = (next_utg->get_front()->is_neighbors_with(exit_node)) ? next_utg->get_front() : next_utg->get_back();
            bool next_reverse = (next_utg->get_front() == next_utg->get_back()) ? single_node_reverse(next_node, exit_node, true) : (next_node == next_utg->get_back());
            next_utg->get_hpc_sequence(walk_seq, next_reverse, true);
            exit_node = next_utg->opposite_end(next_node);
        }
        ret_seqs.emplace_back(walk_seq);
    }
    return true;
}

void phase::build_walk(Walk& walk, Node* start_node, int hap_num){
    int other_hap_num = (hap_num == HAP1) ? HAP2 : HAP1;
    Node* curr = start_node;
    while (curr != NULL){
        Unitig* add_utg = curr->get_unitig();
        int times_walked = add_utg->get_times_walked();
        if (times_walked >= add_utg->get_walk_limit()) break;
        add_utg->set_times_walked(times_walked + 1);
        walk.utgs.push_back(add_utg);
        if (add_utg->get_assignment() == hap_num || add_utg->get_assignment() == HOM){
            walk.num_hap_nodes += add_utg->num_nodes();
        }
        Node* opp_node = add_utg->opposite_end(curr);
        const EdgeList& opp_neighbors = opp_node->view_neighbors();
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

        if (!hap_neighs.empty()){
            choices = &hap_neighs;
        }else if(!other_neighs.empty()){
            choices = &other_neighs;
        }
        if (choices == NULL) break;
        Node* next = (*choices)[0];
        if (choices->size() > 1){
            int best_score {-1};
            for (size_t i{0}; i < choices->size(); i++){
                int budget {LOOKAHEAD_BUDGET};
                int score = phase::lookahead_score((*choices)[i], hap_num, budget);
                if (score > best_score){ best_score = score; next = (*choices)[i]; }
            }
        }
        curr = next;
    }
}

int phase::lookahead_score(Node* to_add, int hap_num, int& budget){
    int other_hap_num = (hap_num == HAP1) ? HAP2 : HAP1;
    Unitig* add_utg = to_add->get_unitig();
    int times_walked = add_utg->get_times_walked();
    if (budget <= 0 || times_walked >= add_utg->get_walk_limit()) return 0;
    budget--;

    add_utg->set_times_walked(times_walked + 1);
    int hap_weight = (add_utg->get_assignment() == hap_num || add_utg->get_assignment() == HOM) ? add_utg->num_nodes() : 0;

    Node* opp_node = add_utg->opposite_end(to_add);
    const EdgeList& opp_neighbors = opp_node->view_neighbors();
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
    std::vector<Node*>* choices = NULL;
    if (!hap_neighs.empty()) {
        choices = &hap_neighs;
    } else if (!other_neighs.empty()) {
        choices = &other_neighs;
    }

    int best_ext {0};
    if (choices != NULL){
        for (size_t i{0}; i < choices->size() && budget > 0; i++){
            int ext = phase::lookahead_score((*choices)[i], hap_num, budget);
            if (ext > best_ext) best_ext = ext;
        }
    }
    add_utg->set_times_walked(times_walked);
    return hap_weight + best_ext;
}

void phase::phase_component(int comp_num, std::vector<std::string>& phased_walks, Graph* prev_graph, std::vector<SketchedSeq*>& corrected_seqs, bool verbose){
    // this will be filled by next function call
    std::vector<Unitig*> utg_list;
    phase::initialize_unitigs_for_tagging(prev_graph, corrected_seqs, utg_list, comp_num);
    if (utg_list.size() == 1){
        std::string utg_seq;
        utg_list[0]->get_hpc_sequence(utg_seq, false, false);
        phased_walks.emplace_back(utg_seq);
    }else{
        bool has_bubble {false};
        for (size_t i{0}; i < utg_list.size(); i++){
            if (utg_list[i]->get_parallel_utg() != NULL){ has_bubble = true; break; }
        }
        if (!has_bubble){
            for (size_t i{0}; i < utg_list.size(); i++) utg_list[i]->set_assignment(HOM);
            std::vector<std::string> path_seqs;
            std::vector<Unitig*> walked;
            if (phase::get_haplotype_sequences(prev_graph, path_seqs, comp_num, &walked) && !path_seqs.empty()){
                phased_walks.emplace_back(path_seqs[0]);
            }
            std::unordered_set<Unitig*> visited(walked.begin(), walked.end());
            size_t n_extra {0};
            for (size_t i{0}; i < utg_list.size(); i++){
                if (visited.count(utg_list[i]) == 0){
                    std::string utg_seq;
                    utg_list[i]->get_hpc_sequence(utg_seq, false, false);
                    phased_walks.emplace_back(utg_seq);
                    n_extra++;
                }
            }
            if (verbose) std::cout << "(" << comp_num << ") no bubble: emitted the homozygous path once plus " << n_extra << " unvisited unitigs" << std::endl;
            return;
        }
        int round {0};
        while(utg_list.size() > 0){
            if (!phase::seed_initial_haplotypes(utg_list)){
                // no more bubbles, remaining unitigs are homozygous
                for (size_t i{0}; i < utg_list.size(); i++){
                    utg_list[i]->set_assignment(HOM);
                }
                break;
            }
            int it {0};
            while(phase::reassign_partitions(prev_graph, corrected_seqs, utg_list) && it < 100){
                if (verbose && it % 10 == 0) std::cout << "(" << comp_num << ") Phasing round " << round << ", iteration " << it << std::endl;
                it++;
            }
            phase::remove_tagged_unitigs(utg_list);
            round++;
        }
        if (verbose) std::cout << "(" << comp_num << ") Phasing assignments complete" << std::endl;
        phase::get_haplotype_sequences(prev_graph, phased_walks, comp_num);
    }
}
