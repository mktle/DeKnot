#include "deknot.h"
#include "graph.h"
#include "syncmer.h"
#include "utils.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>

Graph::Graph(int set_k_prime){
    assert(this->node_table.size() == 0);
    this->k_prime = set_k_prime;
}

Graph::~Graph(){
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        if (this->unitig_list[i] != NULL){
            delete this->unitig_list[i];
        }
    }
    for (size_t i{0}; i < this->node_list.size(); i++){
        if (this->node_list[i] != NULL){
            delete this->node_list[i];
        }
    }
}

int Graph::get_k_prime(){
    return this->k_prime;
}

void Graph::set_k_prime(int set_val){
    // cannot change to a different value if the graph has already started to be constructed
    assert(this->node_table.size() == 0);
    this->k_prime = set_val;
}

int Graph::num_nodes(){
    return (int)this->node_list.size() - (int)this->node_free_list.size();
}

int Graph::num_unitigs(){
    return (int)this->unitig_list.size() - (int)this->unitig_free_list.size();
}

//  run *after* graph is completed (i.e., in the next iteration)
// check if the sequence is completely contained in a unitig from this previous iteration
bool Graph::check_sequence(SketchedSeq* orig_seq, size_t check_k_prime, int& min_span, SketchedSeq*& corrected_seq){
    ankerl::unordered_dense::map<int, int> utg_spans;
    std::vector<DroppedBlock> to_correct;
    size_t check_thresh = orig_seq->kmers.size() - check_k_prime + 1;
    size_t i {0};
    /****************************************
     * Navigate to first non-null k-syncmer *
     ****************************************/

    Span<uint64_t> curr_kmers(orig_seq->kmers.data(), check_k_prime);
    char curr_orient = syncmer::canonical_syncmer(curr_kmers);
    KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};

    // go to the first node that was removed from the graph
    Node* curr_node = this->get_node_from_lookup_table(curr_canon);
    while (curr_node == NULL && i < check_thresh){
        curr_kmers = Span(orig_seq->kmers.data() + i, check_k_prime);
        curr_orient = syncmer::canonical_syncmer(curr_kmers);
        curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};
        curr_node = this->get_node_from_lookup_table(curr_canon);
        i++;
    }
    // where we will start copying sequence from (i.e., drop errors on ends of reads)
    //size_t curr_orig_idx = i == 1 ? 0 : i;
    size_t curr_orig_idx = i;

    bool in_dropped {false};
    Node* dropped_start_bound = NULL;
    Node* prev_node = curr_node;
    int dropped_start_idx {-1};

    /************************************************
     * Identify parts of the read that were dropped *
     ************************************************/
    for(; i < check_thresh; i++){
        Span<uint64_t> kmers(orig_seq->kmers.data() + i, check_k_prime);
        char orient = syncmer::canonical_syncmer(kmers);
        KmerView canon = {kmers, orient == REVERSE_ORIENTATION};
        curr_node = this->get_node_from_lookup_table(canon);

        if (curr_node != NULL){
            if (in_dropped){
                assert(dropped_start_idx != -1);
                // at the end of a dropped block, so add to our list
                if (dropped_start_bound != curr_node){
                    if(dropped_start_bound->get_unitig() == curr_node->get_unitig() || dropped_start_bound->get_unitig()->is_neighbors_with(curr_node->get_unitig())){
                        to_correct.emplace_back(dropped_start_bound, curr_node, dropped_start_idx, i - 1, dropped_start_bound->get_unitig() != curr_node->get_unitig());
                    }
                }
                in_dropped = false;
                dropped_start_idx = -1;
            }
            utg_spans[curr_node->get_unitig()->get_id()]++;
        }else if (!in_dropped)r
            //on a new section of dropped region
            assert(prev_node != NULL);
            in_dropped = true;
            dropped_start_bound = prev_node;
            dropped_start_idx = (int)i;
        }
        prev_node = curr_node;
    }

    // check if we need to drop other end of read because it contains an error
    size_t read_stop_idx = in_dropped ? (size_t)dropped_start_idx : orig_seq->kmers.size();

    /**************************************
     * Correction on underlying sequences *
     **************************************/
    if (to_correct.size() > 0 || in_dropped || curr_orig_idx != 0 || read_stop_idx != orig_seq->kmers.size()){
        corrected_seq = new SketchedSeq(orig_seq->name + "M", true);
        if (orig_seq->name == "m84039_230117_233243_s1/105583716/ccs" || orig_seq->name == "m84039_230117_233243_s1/105583716/ccsM" || orig_seq->name == "m84039_230117_233243_s1/105583716/ccsMM" || orig_seq->name == "m84039_230117_233243_s1/105583716/ccsMMM" || orig_seq->name == "m84039_230117_233243_s1/105583716/ccsMMMM"){//mkl
            std::cout << "CHANGING THIS READ " << orig_seq->name << std::endl;//mkl
        }//mkl
        for (size_t j{0}; j < to_correct.size(); j++){
            // copy over the unaltered sequences
            while (curr_orig_idx <= to_correct[j].start_idx + check_k_prime - 2){
                corrected_seq->kmers.push_back(orig_seq->kmers[curr_orig_idx]);
                corrected_seq->flipped.push_back(orig_seq->flipped[curr_orig_idx]);
                corrected_seq->unique_seq_len.push_back(orig_seq->unique_seq_len[curr_orig_idx]);
                curr_orig_idx++;
            }
            if (!to_correct[j].across_utgs){
                std::cout << "within unitig " << to_correct[j].end1->get_ << std::endl;//mkl
                this->correct_read_error_within_unitig(to_correct[j], corrected_seq);
            }else if (to_correct[j].end1->get_unitig()->is_neighbors_with(to_correct[j].end2->get_unitig())){
                this->correct_read_error_across_unitigs(to_correct[j], corrected_seq);
            }
            curr_orig_idx = to_correct[j].end_idx + check_k_prime;
        }
        while (curr_orig_idx < read_stop_idx){
            corrected_seq->kmers.push_back(orig_seq->kmers[curr_orig_idx]);
            corrected_seq->flipped.push_back(orig_seq->flipped[curr_orig_idx]);
            corrected_seq->unique_seq_len.push_back(orig_seq->unique_seq_len[curr_orig_idx]);
            curr_orig_idx++;
        }
    }else{
        corrected_seq = NULL;
    }
    min_span = -1;
    if (orig_seq->is_read){
        for (auto const& [id, span] : utg_spans) {
            if (min_span == -1 || span <  min_span) min_span = span;
        }
    }
    return utg_spans.size() == 1;
}

void Graph::correct_read_error_within_unitig(DroppedBlock& dropped_block, SketchedSeq*& corrected_seq){
    // now we're in a dropped block, so retrieve the kmers for replacement
    // TODO: linear search is way too slow
    const std::vector<Node*> utg_nodes = dropped_block.end1->get_unitig()->view_nodes();
    int first_idx {-1};
    int last_idx {-1};
    for (size_t k{0}; k < utg_nodes.size(); k++){
        if (first_idx != -1){
            if (utg_nodes[k] == dropped_block.end2) break;
            // within the region we are copying
            bool flipped;
            int unique_sl;
            corrected_seq->kmers.push_back(utg_nodes[k]->unique_kmer_from_neigh(utg_nodes[k - 1], flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
        }else {
            if (utg_nodes[k] == dropped_block.end1){
                first_idx = (int)k;
                // already saw the other bookend node, so terminate
                if (last_idx != -1) break;
            }else if (utg_nodes[k] == dropped_block.end2){
                last_idx = (int)k;
            }
        }
    }
    // if the sequence was backwards, add the relevant part now
    if (last_idx != -1){
        assert(last_idx < first_idx);
        for (int k{first_idx - 1}; k > last_idx; k--){
            bool flipped {false};
            int unique_sl {-1};
            corrected_seq->kmers.push_back(utg_nodes[(size_t)k]->unique_kmer_from_neigh(utg_nodes[(size_t)k + 1], flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
        }
    }
}

void Graph::correct_read_error_across_unitigs(DroppedBlock& dropped_block, SketchedSeq*& corrected_seq){
    Unitig* utg1 = dropped_block.end1->get_unitig();
    Unitig* utg2 = dropped_block.end2->get_unitig();
    assert(utg1 != utg2);

    /********************************************************
     * Determine which parts of each unitig are on the walk *
     ********************************************************/
    bool utg1_frontside;
    bool utg2_frontside;
    
    if (utg1->get_front()->is_neighbors_with(utg2->get_front())){
        utg1_frontside = true;
        utg2_frontside = true;
    }else if(utg1->get_front()->is_neighbors_with(utg2->get_back())){
        utg1_frontside = true;
        utg2_frontside = false;
    }else if(utg1->get_back()->is_neighbors_with(utg2->get_front())){
        utg1_frontside = false;
        utg2_frontside = true;
    }else if(utg1->get_back()->is_neighbors_with(utg2->get_back())){
        utg1_frontside = false;
        utg2_frontside = false;
    }else{
        assert(0);
    }

    /***********************************
     * Add kmers from the first unitig *
     ***********************************/
    const std::vector<Node*> utg1_nodes = utg1->view_nodes();
    int k {0};
    bool flipped {false};
    int unique_sl {-1};
    // find the endpoint in the unitig
    while(utg1_nodes[(size_t)k] != dropped_block.end1){
        k++;
    }
    if (utg1_frontside){
        k--;
        while(k >= 0){
            corrected_seq->kmers.push_back(utg1_nodes[(size_t)k]->unique_kmer_from_neigh(utg1_nodes[(size_t)k+1], flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
            k--;
        }
    }else{
        k++;
        while(k < utg1_nodes.size()){
            corrected_seq->kmers.push_back(utg1_nodes[(size_t)k]->unique_kmer_from_neigh(utg1_nodes[(size_t)k-1], flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
            k++;
        }
    }
    /************************************
     * Add kmers from the second unitig *
     ************************************/
    Node* last_node = utg1_frontside ? utg1->get_front() : utg1->get_back();
    const std::vector<Node*> utg2_nodes = utg2->view_nodes();
    if (utg2_frontside){
        k = 0;
        while(utg2_nodes[(size_t)k]->get_node_id() != dropped_block.end2->get_node_id()){
            corrected_seq->kmers.push_back(utg2_nodes[(size_t)k]->unique_kmer_from_neigh(last_node, flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
            last_node = utg2_nodes[(size_t)k];
            k++;
        }
    }else{
        k = (int)utg2_nodes.size() - 1;
        while(utg2_nodes[(size_t)k]->get_node_id() != dropped_block.end2->get_node_id()){
            corrected_seq->kmers.push_back(utg2_nodes[(size_t)k]->unique_kmer_from_neigh(last_node, flipped, unique_sl));
            corrected_seq->flipped.push_back(flipped);
            corrected_seq->unique_seq_len.push_back(unique_sl);
            last_node = utg2_nodes[(size_t)k];
            k--;
        }
    }
}

bool Graph::add_sequence(SketchedSeq& input_seq, Graph* prev_graph){
    assert(input_seq.kmers.size() >= this->k_prime);
    Span<uint8_t> first_flipped(input_seq.flipped.data(), (size_t)this->k_prime);
    Span<uint8_t> second_flipped;
    Span<int> first_unique_seq_len(input_seq.unique_seq_len.data() + 1, (size_t)this->k_prime - 1);
    Span<int> second_unique_seq_len;
    Span<uint64_t> first_kmers(input_seq.kmers.data(), (size_t)this->k_prime);
    Span<uint64_t> second_kmers;
    char first_orient = syncmer::canonical_syncmer(first_kmers);
    char second_orient;
    KmerView first_ksyncmer = {first_kmers, first_orient == REVERSE_ORIENTATION};
    KmerView second_ksyncmer;

    bool first_connection {true};
    size_t num_ksync_pairs = input_seq.kmers.size() - (size_t)this->k_prime;
    for(size_t i{0}; i < num_ksync_pairs; i++){
        second_flipped = Span(input_seq.flipped.data() + i + 1, (size_t)this->k_prime);
        second_unique_seq_len = Span(input_seq.unique_seq_len.data() + i + 2, (size_t)this->k_prime - 1);
        second_kmers = Span(input_seq.kmers.data() + i + 1, (size_t)this->k_prime);
        second_orient = syncmer::canonical_syncmer(second_kmers);
        KmerView second_ksyncmer = {second_kmers, second_orient == REVERSE_ORIENTATION};

        if (first_connection){
            int increment_val {1};
            if (!input_seq.is_read){
                increment_val = prev_graph->k_plus_one_mer_count(first_kmers);
            }
            this->increment_k_min_mer_count(increment_val, first_ksyncmer, first_flipped, first_unique_seq_len);
            first_connection = false;
        }
        int increment_val {1};
        if (!input_seq.is_read){
            increment_val = prev_graph->k_plus_one_mer_count(second_kmers);
        }
        this->increment_k_min_mer_count(increment_val, second_ksyncmer, second_flipped, second_unique_seq_len);

        // add connection between these adjacent k-min-mers
        this->add_connection(first_ksyncmer, first_orient, second_ksyncmer, second_orient);
        // start shift to next pair
        first_flipped = second_flipped;
        first_unique_seq_len = second_unique_seq_len;
        first_kmers = second_kmers;
        first_orient = second_orient;
        first_ksyncmer = second_ksyncmer;
    }
    return input_seq.is_read;
}


bool Graph::add_connection(KmerView& first, char first_orient, KmerView& second, char second_orient){
    Node* first_node {this->get_node_from_lookup_table(first)};
    Node* second_node {this->get_node_from_lookup_table(second)};
    assert(first_node != NULL);
    assert(second_node != NULL);

    // first check if they're already neighbors—if yes, then everything is already handled
    // TODO: need to actually check the edge direction
    if (first_node->is_neighbors_with(second_node)){
        assert(second_node->is_neighbors_with(first_node));
        return false;
    }

    first_node->add_neighbor(second_node, first_orient, second_orient);
    second_node->add_neighbor(first_node, utils::reverse_orientation(second_orient), utils::reverse_orientation(first_orient));

    // TODO: change so it's not a new unitig every time a node is created
    // will have to make this function handle more cases
    assert(first_node->is_assigned_to_unitig());
    assert(second_node->is_assigned_to_unitig());

    Unitig* first_utg {first_node->get_unitig()};
    Unitig* second_utg {second_node->get_unitig()};

    if (first_node == second_node){
        this->handle_self_connection(first_node);
    }else if (first_utg->node_is_endpoint(first_node) && second_utg->node_is_endpoint(second_node) && first_node->num_edges_in_dir(first_orient) == 1 && second_node->num_edges_in_dir(utils::reverse_orientation(second_orient)) == 1){
        // merge unitigs into one
        if (first_utg != second_utg){
            this->merge_unitigs(first_node, second_node);
        }
    }else if (first_utg->node_is_endpoint(first_node) && second_utg->node_is_endpoint(second_node)){
        this->check_endpoints(first_node, second_node);
    }else{
        // we have to break some unitigs apart, and the nodes will be connected but not on the same unitig
        this->check_endpoints(first_node, second_node);
        first_utg = first_node->get_unitig();
        second_utg = second_node->get_unitig();
        if (!(first_utg->node_is_endpoint(first_node))){
            // same_utg handles special case where the nodes are originally on the same unitig
            bool same_utg = (first_utg->get_id() == second_utg->get_id());
            bool second_orig_end = second_utg->node_is_endpoint(second_node);
            Unitig* new_utg = this->split_unitig(first_utg, first_node);
            if (same_utg && !second_orig_end){
                if (first_utg->contains_node(second_node)){
                    this->split_unitig(first_utg, second_node);
                }else{
                    assert(new_utg->contains_node(second_node));
                    this->split_unitig(new_utg, second_node);
                }
            }
            first_utg = first_node->get_unitig();
            second_utg = second_node->get_unitig();
        }
        if (!second_utg->node_is_endpoint(second_node)){
            this->split_unitig(second_utg, second_node);
        }
    }
    return true;
}

void Graph::handle_self_connection(Node* node){
    Unitig* utg = node->get_unitig();
    Unitig* new_utg = NULL;
    if (utg->num_nodes() > 1){
        // isolate node into its own unitig
        if (utg->get_front() != node && utg->get_back() != node){
            new_utg = this->split_unitig(utg, node);
            // case where node is on the end of the new unitig
            if (new_utg->get_back() == node){
                new_utg->pop_back();
                this->create_unitig_for_node(node);
            }else{
                assert(utg->get_front() == node);
                // case where node is at the beginning of the old unitig
                utg->pop_front();
                this->create_unitig_for_node(node);
            }
        }else{
            //either at the beginning or end of utg
            if (utg->get_front() == node){
                utg->pop_front();
                this->create_unitig_for_node(node);
            }else if(utg->get_back() == node){
                utg->pop_back();
                this->create_unitig_for_node(node);
            }
        }
    }
    assert(node->get_unitig()->num_nodes() == 1);
    utg->get_median_freq();
    if (new_utg != NULL) utg->get_median_freq();
    this->free_unitig_if_empty(utg);
    if (new_utg != NULL) this->free_unitig_if_empty(new_utg);
}

void Graph::check_endpoints(Node* first_node, Node* second_node){
    Unitig* first_utg = first_node->get_unitig();
    Unitig* second_utg = second_node->get_unitig();
    // endpoints may need to be isolated into their own unitigs
    if (first_utg->node_is_endpoint(first_node) && first_node->num_edges_in_dir(first_utg->get_in_dir(first_node)) > 1){
        if (first_utg->get_front() == first_node){
            first_utg->pop_front();
            this->free_unitig_if_empty(first_utg);
        }else{
            assert(first_utg->get_back() == first_node);
            first_utg->pop_back();
            this->free_unitig_if_empty(first_utg);
        }
        this->create_unitig_for_node(first_node);
    }

    if (first_node != second_node && second_utg->node_is_endpoint(second_node) && second_node->num_edges_in_dir(second_utg->get_in_dir(second_node)) > 1){
        if (second_utg->get_front() == second_node){
            second_utg->pop_front();
            this->free_unitig_if_empty(second_utg);
        }else{
            assert(second_utg->get_back() == second_node);
            second_utg->pop_back();
            this->free_unitig_if_empty(second_utg);
        }
        this->create_unitig_for_node(second_node);
    }
}

// the nodes of the second unitig will be dumped into the first unitig
void Graph::merge_unitigs(Node* first_node, Node* second_node){
    Unitig* first_utg {first_node->get_unitig()};
    Unitig* second_utg {second_node->get_unitig()};
    assert(first_utg->node_is_endpoint(first_node));
    assert(second_utg->node_is_endpoint(second_node));
    assert(first_utg != second_utg);

    const std::vector<Node*>& second_nodes {second_utg->view_nodes()};
    if (first_utg->get_front()->get_node_id() == first_node->get_node_id()){
        if (second_utg->get_front()->get_node_id() == second_node->get_node_id()){
            for (int i{0}; i < second_utg->num_nodes(); i++){
                first_utg->add_node_front(second_nodes[(size_t)i]);
                second_nodes[(size_t)i]->assign_unitig(first_utg);
            }
        }else{
            assert(second_utg->get_back()->get_node_id() == second_node->get_node_id());
            for (int i{second_utg->num_nodes() - 1}; i >= 0; i--){
                first_utg->add_node_front(second_nodes[(size_t)i]);
                second_nodes[(size_t)i]->assign_unitig(first_utg);
            }
        }
    }else{
        assert(first_utg->get_back()->get_node_id() == first_node->get_node_id());
        if (second_utg->get_front()->get_node_id() == second_node->get_node_id()){
            for (int i{0}; i < second_utg->num_nodes(); i++){
                first_utg->add_node_back(second_nodes[(size_t)i]);
                second_nodes[(size_t)i]->assign_unitig(first_utg);
            }
        }else{
            assert(second_utg->get_back()->get_node_id() == second_node->get_node_id());
            for (int i{second_utg->num_nodes() - 1}; i >= 0; i--){
                first_utg->add_node_back(second_nodes[(size_t)i]);
                second_nodes[(size_t)i]->assign_unitig(first_utg);
            }
        }
    }
    second_utg->clear_nodes();
    this->destroy_unitig(second_utg);
}

Unitig* Graph::split_unitig(Unitig* to_split, Node* breakpoint){
    std::vector<Node*> drained_nodes;
    int new_id {this->next_unitig_id()};

    // create unitig for nodes to the left of the breakpoint
    to_split->drain_nodes(drained_nodes, breakpoint);
    assert(to_split->num_nodes() > 0);
    Unitig* new_utg = new Unitig(drained_nodes, new_id);
    this->add_unitig_to_list(new_utg);
    for (size_t i{0}; i < drained_nodes.size(); i++){
        drained_nodes[i]->assign_unitig(new_utg);
    }
    return new_utg;
}

int Graph::next_node_id(){
    if (this->node_free_list.size() > 0){
        int next{this->node_free_list.back()};
        this->node_free_list.pop_back();
        return next;
    }
    return (int)this->node_list.size();
}

int Graph::next_unitig_id(){
    if (this->unitig_free_list.size() > 0){
        int next{this->unitig_free_list.back()};
        this->unitig_free_list.pop_back();
        return next;
    }
    return (int)this->unitig_list.size();
}

// k_min_mer should be reversed already for lookup
// but the bookkeping structures won't be (to reduce number of copies)
void Graph::increment_k_min_mer_count(int increment, KmerView& k_syncmer, Span<uint8_t>& kmer_flipped, Span<int>& unique_seq_idx){
    Node* node_ptr = get_node_from_lookup_table(k_syncmer);
    if (node_ptr == NULL){
        // have not encountered this k_min_mer before, so need to create a new node and unitig
        Node* new_node = new Node(k_syncmer, kmer_flipped, unique_seq_idx, this->next_node_id(), NULL);
        this->add_node_to_list(new_node);
        this->create_unitig_for_node(new_node);

        this->add_node_to_lookup_table(new_node);
        new_node->set_count(increment);
    }else{
        // have already encountered this k_min_mer, so we just need to increment count
        node_ptr->set_count(node_ptr->get_count() + increment);
    }
}

void Graph::increment_unitig(int increment, Unitig* to_increment){
    const std::vector<Node*> nodes = to_increment->view_nodes();
    for (size_t i{0}; i < nodes.size(); i++){
        nodes[i]->set_count(nodes[i]->get_count() + increment);
    }
}

int Graph::k_plus_one_mer_count(Span<uint64_t>& k_syncmer){
    size_t curr_k_prime {k_syncmer.size()};
    assert(curr_k_prime > this->k_prime);
    // find the minimum value across all sub-(k'-n)-syncmers
    Span<uint64_t> copy = Span(k_syncmer.data(), (size_t)this->k_prime);
    KmerView canon = {copy, syncmer::canonical_syncmer(copy) == REVERSE_ORIENTATION};
    int min_val = this->node_table.count(canon) > 0 ? this->node_table[canon]->get_count() : 1;

    // check all the other (k'-n) syncmers
    for (int i{0}; i < curr_k_prime - (size_t)this->k_prime; i++){
        copy = Span(k_syncmer.data() + i + 1, (size_t)this->k_prime);
        canon = {copy, syncmer::canonical_syncmer(copy) == REVERSE_ORIENTATION};
        int curr_val = this->node_table.count(canon) > 0 ? this->node_table[canon]->get_count() : 1;
        if (curr_val < min_val){
            min_val = curr_val;
        }
    }
    return min_val;
}

void Graph::add_node_to_list(Node* new_node){
    if (new_node->get_node_id() < this->node_list.size()){
        assert(this->node_list[(size_t)new_node->get_node_id()] == NULL);
        this->node_list[(size_t)new_node->get_node_id()] = new_node;
    }else{
        this->node_list.push_back(new_node);
    }
}

void Graph::add_unitig_to_list(Unitig* new_utg){
    if (new_utg->get_id() < this->unitig_list.size()){
        assert(this->unitig_list[(size_t)new_utg->get_id()] == NULL);
        this->unitig_list[(size_t)new_utg->get_id()] = new_utg;
    }else{
        this->unitig_list.push_back(new_utg);
    }
}

void Graph::destroy_unitig(Unitig* to_remove){
    assert(this->unitig_list[(size_t)to_remove->get_id()] != NULL);
    this->unitig_list[(size_t)to_remove->get_id()] = NULL;
    this->unitig_free_list.push_back(to_remove->get_id());
    // free all nodes

    const std::vector<Node*> members = to_remove->view_nodes();
    for (size_t i{0}; i < members.size(); i++){
        this->destroy_node(members[i]);
    }
    delete to_remove;
}

void Graph::destroy_node(Node* to_remove){
    // we should only be trying to delete disconnected or empty unitigs
    // so we shouldn't have to worry about the endpoints still being connected to extant nodes
    this->node_free_list.push_back(to_remove->get_node_id());
    this->node_list[(size_t)to_remove->get_node_id()] = NULL;
    this->node_table.erase(to_remove->get_kmer_view());
    delete to_remove;
}

void Graph::create_unitig_for_node(Node* node){
    int utg_id {this->next_unitig_id()};
    std::vector<Node*> node_vec;
    node_vec.push_back(node);
    Unitig* new_utg = new Unitig(node_vec, utg_id);
    this->add_unitig_to_list(new_utg);
    node->assign_unitig(new_utg);
}

void Graph::trim_tips(){
    bool changed = true;
    while(changed){
        changed = false;
        std::vector<Unitig*> sorted_unitigs;
        for (size_t i{0}; i < this->unitig_list.size(); i++){
            Unitig* utg = this->unitig_list[i];
            if (utg != NULL){
                utg->update_median_freq();
                sorted_unitigs.push_back(utg);
            }
        }
        std::sort(sorted_unitigs.begin(), sorted_unitigs.end(), compare_unitigs);

        while(sorted_unitigs.size() > 0){
            Unitig* utg = sorted_unitigs.back();
            Node* keep_node;
            Node* inner_end = utg->get_tip_inner_end(keep_node);
            if (inner_end != NULL){
                assert(keep_node != NULL);
                Unitig* keep_utg = keep_node->get_unitig();
                if (utg->get_median_freq() < 0.5 * keep_utg->get_median_freq() && utg->num_nodes() < 500){
                    changed = true;
                    this->disconnect_unitig_end(inner_end, keep_utg, sorted_unitigs);
                    this->destroy_unitig(utg);
                }
            }
            sorted_unitigs.pop_back();
        }
    }
}

void Graph::remove_floaters(int threshold){
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* utg = this->unitig_list[i];
        if (utg != NULL && utg->is_floating() && utg->get_median_freq() <= threshold){
            this->destroy_unitig(utg);
        }
    }
}

bool Graph::trim_low_freq_unitigs(int threshold, bool check_dinuc){
    // we will start by merging from longer unitigs
    // so first sort the unitigs by length

    // make sure the median frequency is updated for each unitig
    std::vector<Unitig*> sorted_unitigs;
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* utg = this->unitig_list[i];
        if (utg != NULL){
            utg->update_median_freq();
            sorted_unitigs.push_back(utg);
        }
    }
    std::sort(sorted_unitigs.begin(), sorted_unitigs.end(), compare_unitigs);
    bool graph_changed {false};
    while(sorted_unitigs.size() > 0){
        Unitig* curr_utg = sorted_unitigs.front();
        assert(curr_utg != NULL);
        assert(curr_utg->num_nodes() > 0);

        float current_thresh {2};
        while (current_thresh <= threshold){
            bool compacted_front = this->trim_low_freq_end(curr_utg->get_front(), current_thresh, sorted_unitigs, check_dinuc);
            bool compacted_back;
            if (curr_utg->get_front() == curr_utg->get_back()){
                // only need to handle one side for single node unitigs
                compacted_back = compacted_front;
            }else{
                compacted_back = this->trim_low_freq_end(curr_utg->get_back(), current_thresh, sorted_unitigs, check_dinuc);
            }

            if (compacted_front || compacted_back){
                // if either of the ends extended the unitig after dropping unitigs, then start the compaction process over
                current_thresh = 2;
                graph_changed = true;
            }else{
                current_thresh++;
            }
        }
        sorted_unitigs.erase(sorted_unitigs.begin());
    }
    return graph_changed;
}

bool Graph::trim_low_freq_end(Node* node_end, float thresh, std::vector<Unitig*>& sorted_list, bool check_dinuc){
    Unitig* curr_utg = node_end->get_unitig();
    const std::vector<Edge> neighbors = node_end->view_neighbors();
    float lowest_freq {-1};
    Unitig* lower_utg {NULL};
    int num_checked{0};
    for (size_t i{0}; i < neighbors.size(); i++){
        Node* neigh_node = neighbors[i].neighbor;
        Unitig* neigh_utg = neigh_node->get_unitig();
        if (neigh_utg != curr_utg){
            if (neigh_utg->get_median_freq() < thresh){
                this->drop_unitig(sorted_list, node_end, neigh_utg);
                return true;
            }else{
                num_checked++;
                if (lower_utg == NULL || neigh_utg->get_median_freq() < lowest_freq){
                    lowest_freq = neigh_utg->get_median_freq();
                    lower_utg = neigh_utg;
                }
            }
        }
    }
    // compacted no neighbors based on frequency, so check for branching dinucleotide site
    // and make sure it isn't a main bubble (true variation at a dinuc site)
    if (check_dinuc && node_end->end_is_dinucleotide() && num_checked == 2 && this->dinuc_topology(node_end, lower_utg)){
        std::cout << "Dropping based on dinucleotide criteria (utg " << lower_utg->get_id() << " dropped during check on utg" << curr_utg->get_id() << ")" << std::endl;
        assert(lower_utg != NULL);
        bool compacted_graph = this->drop_unitig(sorted_list, node_end, lower_utg);
        return compacted_graph;
    }
    return false;
}

bool Graph::dinuc_topology(Node* curr_end, Unitig* dropped_utg){
    // ask ourselves if we can get to the dropped_utg's other neighbor without going through the curr_utg's neighbors
    // if so, this is likely a dinucleotide error and not a real variant
    Unitig* curr_utg = curr_end->get_unitig();
    const std::vector<Edge> neighbors = curr_end->view_neighbors();
    ankerl::unordered_dense::set<Node*> seen_nodes;
    seen_nodes.insert(curr_end);
    for (size_t i{0}; i < neighbors.size(); i++){
        Node* neigh_node = neighbors[i].neighbor;
        Unitig* neigh_utg = neigh_node->get_unitig();
        if (neigh_utg != dropped_utg && neigh_utg != curr_utg){
            seen_nodes.insert(neigh_node);
            seen_nodes.insert(neigh_utg->opposite_end(neigh_node));
        }
    }
    // perform a graph search
    std::queue<Node*> to_check;
    to_check.push(curr_utg->opposite_end(curr_end));
    seen_nodes.insert(curr_end);
    while(to_check.size() > 0){
        Node* curr_node {to_check.front()};
        to_check.pop();

        const std::vector<Edge> curr_neighbors = curr_node->view_neighbors();
        for (size_t i{0}; i < curr_neighbors.size(); i++){
            Node* neigh_node = curr_neighbors[i].neighbor;
            Unitig* neigh_utg = neigh_node->get_unitig();
            if (curr_node->get_unitig() != neigh_utg){
                if (neigh_utg == dropped_utg){
                    return true;
                }
                if (seen_nodes.find(neigh_node) == seen_nodes.end()){
                    seen_nodes.insert(neigh_node);
                    to_check.push(neigh_node);
                }
            }
        }
        Node* curr_op = curr_node->get_unitig()->opposite_end(curr_node);
        if (seen_nodes.find(curr_op) == seen_nodes.end()){
            seen_nodes.insert(curr_op);
            to_check.push(curr_op);
        }
    }
    return false;
}

bool Graph::drop_unitig(std::vector<Unitig*>& sorted_list, Node* curr_node, Unitig* neigh_utg){
    Unitig* curr_utg = curr_node->get_unitig();
    auto it = std::find(sorted_list.begin(), sorted_list.end(), neigh_utg);
    if(it != sorted_list.end()){
        sorted_list.erase(it);
    }
    // first handle special case where the dropped unitig is connected to curr unitig on both ends
    if (neigh_utg->num_nodes() > 1 && curr_node->is_neighbors_with(neigh_utg->get_front()) && curr_node->is_neighbors_with(neigh_utg->get_back())){
        curr_node->remove_neighbor(neigh_utg->get_front());
        curr_node->remove_neighbor(neigh_utg->get_back());
        this->destroy_unitig(neigh_utg);
        return false;
    }

    bool compacted_front = this->disconnect_unitig_end(neigh_utg->get_front(), curr_utg, sorted_list);
    bool compacted_back;
    if (neigh_utg->get_front() == neigh_utg->get_back()){
        compacted_back = compacted_front;
    }else{
        compacted_back = this->disconnect_unitig_end(neigh_utg->get_back(), curr_utg, sorted_list);
    }
    this->destroy_unitig(neigh_utg);
    if (compacted_front || compacted_back){
        return true;
    }
    return false;
}

// return true if base_utg gets merged with another unitig
bool Graph::disconnect_unitig_end(Node* remove_node, Unitig* keep_utg, std::vector<Unitig*>& sorted_list){
    Unitig* remove_utg = remove_node->get_unitig();
    assert(remove_utg->node_is_endpoint(remove_node));

    const std::vector<Edge> remove_neighbors = remove_node->view_neighbors();
    // remove the endpoint node from its neighbors' edges
    bool compacted {false};
    for (size_t i{0}; i < remove_neighbors.size(); i++){
        Node* neigh_node = remove_neighbors[i].neighbor;
        // don't need to handle neighbors on the same unitig
        if (neigh_node->get_unitig() != remove_utg){
            Unitig* neigh_utg = neigh_node->get_unitig();
            assert(neigh_utg->node_is_endpoint(neigh_node));

            neigh_node->remove_neighbor(remove_node);

            // since we are in the point of view of the removed node, we have to flip the direction to get the direction from the point of view of the neighbor
            char dropped_dir {utils::reverse_orientation(remove_neighbors[i].neighbor_orientation)};

            const std::vector<Edge> neigh_neighbors = neigh_node->view_neighbors();
            // check if we can now recompact the neighbor unitig with one of its own neighbors
            if (neigh_node->num_edges_in_dir(dropped_dir) == 1){
                // get the neighbor's neighbor with the right edge
                for (size_t j{0}; j < neigh_neighbors.size(); j++){
                    if (neigh_neighbors[j].self_orientation == dropped_dir && neigh_neighbors[j].neighbor->get_unitig() != neigh_utg && neigh_neighbors[j].neighbor->num_edges_in_dir(utils::reverse_orientation(neigh_neighbors[j].neighbor_orientation)) == 1){
                        // need to ensure we don't discard the unitig from the highest level loop (in trim_low_freq_unitigs)
                        Node* donor;
                        Node* recipient;
                        if (neigh_neighbors[j].neighbor->get_unitig() == keep_utg){
                            recipient = neigh_neighbors[j].neighbor;
                            donor = neigh_node;
                        }else{
                            donor = neigh_neighbors[j].neighbor;
                            recipient = neigh_node;
                        }
                        auto it = std::find(sorted_list.begin(), sorted_list.end(), donor->get_unitig());
                        if (it != sorted_list.end()){
                            sorted_list.erase(it);
                        }

                        this->merge_unitigs(recipient, donor);
                        recipient->get_unitig()->update_median_freq();
                        compacted = true;
                    }
                }
            }
        }
    }
    return compacted;
}

void Graph::add_node_to_lookup_table(Node* to_add){
    this->node_table[to_add->get_kmer_view()] = to_add;
}

Node* Graph::get_node_from_lookup_table(KmerView& k_syncmer){
    if (this->node_table.count(k_syncmer) > 0){
        return this->node_table[k_syncmer];
    }
    return NULL;
}

bool Graph::any_connected_unitigs(){
    // check if any unitigs are connected to another unitig
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        if (this->unitig_list[i] != NULL && this->unitig_list[i]->neighbors_unitig()){
            return true;
        }
    }
    return false;
}

void Graph::write_nodes_to_gfa(std::string& file_name, bool verbose){
    std::ofstream out_file(file_name);
    for (size_t i{0}; i < this->node_list.size(); i++){
        if (node_list[i] != NULL){
            out_file << "S\tn" << i << "\t";
            if (verbose){
                std::string sequence {""};
                node_list[i]->get_hpc_sequence(sequence, false);
                out_file << sequence;
            }else{
                out_file << "A";
            }
            out_file << "\tMF:i:" << this->node_list[i]->get_count() << '\n';

            const std::vector<Edge> neighs = this->node_list[i]->view_neighbors();
            for (size_t j{0}; j < neighs.size(); j++){
                out_file << "L\tn" << i << "\t" << neighs[j].self_orientation << "\tn" << neighs[j].neighbor->get_node_id() << "\t" << neighs[j].neighbor_orientation << "\t0M\n";
            }
        }
    }
    out_file.close();
}

void Graph::write_hpc_unitigs_to_gfa(std::string& file_name, bool verbose){
    std::vector<std::string> utg_seqs;
    this->get_all_hpc_seqs(utg_seqs);

    std::ofstream out_file(file_name);
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* curr_utg = this->unitig_list[i];
        if (curr_utg != NULL){
            assert(curr_utg->num_nodes() > 0);
            curr_utg->update_median_freq();
            out_file << "S\tutg" << curr_utg->get_id() << '\t';

            if (verbose){
                out_file << utg_seqs[i];
            }else{
                out_file << "*\tLN:i:" << utg_seqs[i].size();
            }
            out_file << "\tMF:i:" << curr_utg->get_median_freq() << '\n';

            if (curr_utg->num_nodes() > 1){
                Node* front_node = curr_utg->get_front();
                const std::vector<Edge> front_neighs = front_node->view_neighbors();
                for (size_t m{0}; m < front_neighs.size(); m++){
                    Node* neigh_node = front_neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_id() != curr_utg->get_id()){
                        char write_dir;
                        if (neigh_utg->num_nodes() == 1){
                            write_dir = front_neighs[m].neighbor_orientation;
                        }else if (neigh_utg->get_front() == neigh_node){
                            write_dir = FORWARD_ORIENTATION;
                        }else{
                            assert(neigh_utg->get_back() == neigh_node);
                            write_dir = REVERSE_ORIENTATION;
                        }
                        out_file << "L\tutg" << curr_utg->get_id() << "\t-\tutg" << neigh_utg->get_id() << "\t" << write_dir << "\t0M\n";
                    }
                }

                Node* back_node = curr_utg->get_back();
                const std::vector<Edge> back_neighs = back_node->view_neighbors();
                for (size_t m{0}; m < back_neighs.size(); m++){
                    Node* neigh_node = back_neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_id() != curr_utg->get_id()){
                        char write_dir;
                        if (neigh_utg->num_nodes() == 1){
                            write_dir = back_neighs[m].neighbor_orientation;
                        }else if (neigh_utg->get_front() == neigh_node){
                            write_dir = FORWARD_ORIENTATION;
                        }else{
                            assert(neigh_utg->get_back() == neigh_node);
                            write_dir = REVERSE_ORIENTATION;
                        }
                        out_file << "L\tutg" << curr_utg->get_id() << "\t+\tutg" << neigh_utg->get_id() << "\t" << write_dir << "\t0M\n";
                    }
                }
            }else{
                // if only one node in unitig, assign in and out edges using orientation
                const std::vector<Edge> neighs = curr_utg->get_front()->view_neighbors();
                for (size_t m{0}; m < neighs.size(); m++){
                    Node* neigh_node = neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_front() == neigh_node){
                        out_file << "L\tutg" << curr_utg->get_id() << "\t" << neighs[m].self_orientation << "\tutg" << neigh_utg->get_id() << "\t+\t0M\n";
                    }else{
                        assert(neigh_utg->get_back() == neigh_node);
                        out_file << "L\tutg" << curr_utg->get_id() << "\t" << neighs[m].self_orientation << "\tutg" << neigh_utg->get_id() << "\t-\t0M\n";
                    }
                }
            }
            out_file << "# n" << curr_utg->get_front()->get_node_id() << "\tn" << curr_utg->get_back()->get_node_id() << '\n';
        }
    }
    out_file.close();
}

void Graph::free_unitig_if_empty(Unitig* to_check){
    if (to_check->num_nodes() == 0){
        this->destroy_unitig(to_check);
    }
}

void Graph::print_nodes(bool verbose){
    for (size_t i{0}; i < this->node_list.size(); i++){
        if (node_list[i] != NULL){
            std::cout << "(" << this->node_list[i]->get_node_id() << ")[" << this->node_list[i]->get_unitig()->get_id() << "]{" << &(this->node_list[i]) << "} " ;
            const std::vector<Edge> neighs = this->node_list[i]->view_neighbors();
            for (size_t j{0}; j < neighs.size(); j++){
                std::cout << neighs[j].neighbor->get_node_id() << ' ';
            }
            std::cout << std::endl;
            if (verbose){
                std::string sequence {""};
                this->node_list[i]->get_hpc_sequence(sequence, false);
                std::cout << sequence << std::endl;
                sequence = "";
                this->node_list[i]->get_hpc_sequence(sequence, true);
                std::cout << sequence << std::endl;
            }
        }
    }
}

void Graph::print_unitigs(bool verbose, bool print_seq){
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        if (unitig_list[i] != NULL){
            if (verbose){
                unitig_list[i]->print(false);
                if (print_seq){
                    std::string sequence {""};
                    this->unitig_list[i]->get_hpc_sequence(sequence, false, false);
                    std::cout << sequence << std::endl;
                }
            }else{
                std::cout << unitig_list[i]->get_id() << ' ';
            }
        }
    }
    if (!verbose){
        std::cout << std::endl;
    }
}

void Graph::get_all_hpc_seqs(std::vector<std::string>& hpc_seqs){
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* utg = this->unitig_list[i];
        std::string sequence{""};
        if (utg != NULL){
            utg->get_hpc_sequence(sequence, false, false);
        }
        // need to push back an empty string for deleted unitigs as well so indices line up later
        hpc_seqs.push_back(std::move(sequence));
    }
}

void Graph::add_unitigs_to_input(std::vector<SketchedSeq*>& input_seqs, int s){
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* utg = this->unitig_list[i];
        if (utg != NULL){
            std::string seq_str {""};
            utg->get_hpc_sequence(seq_str, false, false);
            std::vector<int> indices;
            std::string name {"utg" + std::to_string(i)};
            SketchedSeq* curr_seq = new SketchedSeq(name, false);
            input_seqs.push_back(curr_seq);
            syncmer::get_all_syncmers(seq_str, s, curr_seq->kmers, curr_seq->flipped, curr_seq->unique_seq_len, indices);
        }
    }
}

const std::vector<Unitig*>& Graph::view_unitigs(){
    return this->unitig_list;
}

dk_graph_t* Graph::convert_to_gfa_representation(){
    size_t num_null {0};
    // need to remap IDs because of null entries in unitig list
    ankerl::unordered_dense::map<int, int> new_utg_ids;
    dk_graph_t* gfa_rep = new dk_graph_t();
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        if (this->unitig_list[i] == NULL){
            num_null++;
        }else{
            new_utg_ids[(int)i] = (int)(i - num_null);
        }
    }

    std::vector<dk_edge_t> edge_list;
    gfa_rep->unitigs = new char*[this->unitig_list.size() - num_null];
    gfa_rep->n_unitigs = (int)(this->unitig_list.size() - num_null);

    /************************************************************
     * Copy over unitig sequences and calculate number of edges *
     ************************************************************/
    for (size_t i{0}; i < this->unitig_list.size(); i++){
        Unitig* curr_utg = this->unitig_list[i];
        if (curr_utg != NULL){
            assert(curr_utg->num_nodes() > 0);
            // copy over unitig sequence 
            int new_id = new_utg_ids[(int)i];
            std::string sequence;
            curr_utg->get_hpc_sequence(sequence, false, false);
            gfa_rep->unitigs[new_id] = new char[sequence.size() + 1];
            std::strcpy(gfa_rep->unitigs[new_id], sequence.c_str());

            // now find edges
            curr_utg->update_median_freq();
            if (curr_utg->num_nodes() > 1){
                Node* front_node = curr_utg->get_front();
                const std::vector<Edge> front_neighs = front_node->view_neighbors();
                for (size_t m{0}; m < front_neighs.size(); m++){
                    Node* neigh_node = front_neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_id() != curr_utg->get_id()){
                        char write_dir;
                        if (neigh_utg->num_nodes() == 1){
                            write_dir = front_neighs[m].neighbor_orientation;
                        }else if (neigh_utg->get_front() == neigh_node){
                            write_dir = FORWARD_ORIENTATION;
                        }else{
                            assert(neigh_utg->get_back() == neigh_node);
                            write_dir = REVERSE_ORIENTATION;
                        }
                        edge_list.push_back({new_id, new_utg_ids[neigh_utg->get_id()], '-', write_dir});
                    }
                }

                Node* back_node = curr_utg->get_back();
                const std::vector<Edge> back_neighs = back_node->view_neighbors();
                for (size_t m{0}; m < back_neighs.size(); m++){
                    Node* neigh_node = back_neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_id() != curr_utg->get_id()){
                        char write_dir;
                        if (neigh_utg->num_nodes() == 1){
                            write_dir = back_neighs[m].neighbor_orientation;
                        }else if (neigh_utg->get_front() == neigh_node){
                            write_dir = FORWARD_ORIENTATION;
                        }else{
                            assert(neigh_utg->get_back() == neigh_node);
                            write_dir = REVERSE_ORIENTATION;
                        }
                        edge_list.push_back({new_id, new_utg_ids[neigh_utg->get_id()], '+', write_dir});
                    }
                }
            }else{
                // if only one node in unitig, assign in and out edges using orientation
                const std::vector<Edge> neighs = curr_utg->get_front()->view_neighbors();
                for (size_t m{0}; m < neighs.size(); m++){
                    Node* neigh_node = neighs[m].neighbor;
                    Unitig* neigh_utg = neigh_node->get_unitig();
                    if (neigh_utg->get_front() == neigh_node){
                        edge_list.push_back({new_utg_ids[curr_utg->get_id()], new_utg_ids[neigh_utg->get_id()], neighs[m].self_orientation, '+'});
                    }else{
                        assert(neigh_utg->get_back() == neigh_node);
                        edge_list.push_back({new_utg_ids[curr_utg->get_id()], new_utg_ids[neigh_utg->get_id()], neighs[m].self_orientation, '-'});
                    }
                }
            }
        }
    }
    gfa_rep->edges = new dk_edge_t[edge_list.size()];
    std::memcpy(gfa_rep->edges, edge_list.data(), edge_list.size() * sizeof(dk_edge_t));
    gfa_rep->n_edges = (int)edge_list.size();
    return gfa_rep;
}
